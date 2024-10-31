import os
import json
import pandas as pd
import sys
import argparse

def process_star_amr(star_amr_path):
    results = {
        "mlst_result": {},
        "plasmid_finder": [],
        "resfinder_genes": []
    }
    tsv_file = os.path.join(star_amr_path, "detailed_summary.tsv")
    #print (tsv_file)
    #sys.exit()
    if os.path.exists(tsv_file):
        df = pd.read_csv(tsv_file, sep="\t")
        for _, row in df.iterrows():
            data_type = row["Data Type"]
            data = row["Data"]
            if data_type == "MLST":
                mlst, scheme = data.split(" ")
                results["mlst_result"]["mlst_sequence"] = mlst
                results["mlst_result"]["mlst_scheme"] = scheme.strip("()")
            elif data_type == "Plasmid":
                results["plasmid_finder"].append(data)
            elif data_type == "Resistance":
                gene = data
                phenotypes = row["Predicted Phenotype"].split(", ")
                results["resfinder_genes"].append({"gene": gene, "phenotypes": phenotypes})
    return results

def process_abricate(abricate_path):
    vfdb_results = []
    tsv_file = os.path.join(abricate_path, "amr.vfdb.results.tsv")
    if os.path.exists(tsv_file):
        df = pd.read_csv(tsv_file, sep="\t")
        for _, row in df.iterrows():
            gene = row["GENE"]
            product_resistance = row["PRODUCT"]
            vfdb_results.append({"gene": gene, "product_resistance": product_resistance})
    return vfdb_results

def process_merge(merge_path):
    mob_rgi_results = []
    csv_file = os.path.join(merge_path, "merged_tables.csv")
    print(csv_file)
    sys.exit()
    if os.path.exists(csv_file):
        df = pd.read_csv(csv_file)
        for _, row in df.iterrows():
            result = {
                "cut_off": row["Cut_Off"],
                "best_hit_aro": row["Best_Hit_ARO"],
                "model_type": row["Model_type"],
		"drug_class": row["Drug Class"].split(";"),
                "resistance_mechanism": row["Resistance Mechanism"].split(";"),
                "amr_gene_families": row["AMR Gene Family"].split(";"),
                "mob_suite_results": {
                    "molecule_type": row["molecule_type"],
                    "primary_cluster_id": row["primary_cluster_id"],
                    "secondary_cluster_id": row["secondary_cluster_id"],
                    "amr_relaxase_type": row["relaxase_type(s)"].split(","),
                    "amr_mpf_type": row["mpf_type"].split(","),
                    "amr_orit_type": row["orit_type(s)"].split(","),
                    "amr_predicted_mobility": row["predicted_mobility"].split(",")
                }
            }
            mob_rgi_results.append(result)
    return mob_rgi_results

def process_ectyper(ectyper_path):
    ectyper_data = {}
    tsv_file = os.path.join(ectyper_path, "output.tsv")
    if os.path.exists(tsv_file):
        df = pd.read_csv(tsv_file, sep="\t")
        ectyper_data = {
            "serotype": df["Serotype"].iloc[0],
            "htype": df["H-type"].iloc[0],
            "otype": df["O-type"].iloc[0]
        }
    return ectyper_data

def process_sistr(sistr_path):
    sistr_data = {}
    tsv_file = os.path.join(sistr_path, "sistr.tab")
    if os.path.exists(tsv_file):
        df = pd.read_csv(tsv_file, sep="\t")
        sistr_data = {
            "serovar": df["serovar"].iloc[0],
            "serogroup": df["serogroup"].iloc[0],
            "serovar_antigen": df["serovar_antigen"].iloc[0],
            "o_antigen": df["o_antigen"].iloc[0],
            "h1": df["h1"].iloc[0],
            "h2": df["h2"].iloc[0],
            "genome": df["genome"].iloc[0]
        }
    return sistr_data

def process_virulencefinder(virulencefinder_path):
    virulence_vf = []
    json_file = os.path.join(virulencefinder_path, "data.json")
    if os.path.exists(json_file):
        with open(json_file) as f:
            data = json.load(f)
            #print(data)
            #sys.exit()
            for species, dbs in data["virulencefinder"]["results"].items():
                for db_name, contigs in dbs.items():
                    for contig, result in contigs.items():
                        vf_gene = result.get("virulence_gene")
                        vf_protein_function = result.get("protein_function")
                        if vf_gene and vf_protein_function:
                            virulence_vf.append({"vf_gene": vf_gene, "vf_protein_function": vf_protein_function})
    return virulence_vf

def main(results_dir):
    results = {}
    
    for isolate_id in os.listdir(results_dir):
        isolate_dir = os.path.join(results_dir, isolate_id)
        if not os.path.isdir(isolate_dir):
            continue
        
        results[isolate_id] = {}
        
        # Process each subdirectory if it exists
        staramr_path = os.path.join(isolate_dir, "STARAMR", "out")
        #print(staramr_path)
        #sys.exit()
        abricate_path = os.path.join(isolate_dir, "ABRICATE")
        merge_path = os.path.join(isolate_dir, "Merge")
        ectyper_path = os.path.join(isolate_dir, "ECTYPER", "out")
        sistr_path = os.path.join(isolate_dir, "SISTR")
        virulencefinder_path = os.path.join(isolate_dir, "VIRULENCEFINDER", "out")
        
        if os.path.exists(staramr_path):
           
            results[isolate_id]["staramr"] = process_star_amr(staramr_path)
        
        if os.path.exists(abricate_path):
            results[isolate_id]["vfdb_results"] = process_abricate(abricate_path)
        
        if os.path.exists(merge_path):
            results[isolate_id]["mob_rgi_results"] = process_merge(merge_path)
        
        if os.path.exists(ectyper_path):
            results[isolate_id]["ectyper"] = process_ectyper(ectyper_path)
        
        if os.path.exists(sistr_path):
            results[isolate_id]["sistr"] = process_sistr(sistr_path)
        
        if os.path.exists(virulencefinder_path):
            results[isolate_id]["virulencefinder"] = process_virulencefinder(virulencefinder_path)
    
    # Save results to JSON
    output_file = "results_summary.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=4)

    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process multiple isolates in the results directory.")
    parser.add_argument("results_dir", type=str, help="Path to the results directory containing isolate subdirectories.")
    
    args = parser.parse_args()
    main(args.results_dir)


