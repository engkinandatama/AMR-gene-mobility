import os
import glob
import pandas as pd

def main():
    bench_dir = r"E:\amien_\Downloads\AMR\pilot_run\benchmarks"
    if not os.path.exists(bench_dir):
        print(f"Benchmark directory not found: {bench_dir}")
        return
    
    files = glob.glob(os.path.join(bench_dir, "*.tsv"))
    records = []
    
    for f in files:
        basename = os.path.basename(f)
        
        # Determine the rule name based on the prefix of the file name
        rule_name = None
        if basename.startswith("host_removal_"):
            rule_name = "host_removal"
        elif basename.startswith("aggregate_"):
            rule_name = "aggregate"
        elif basename.startswith("download_"):
            rule_name = "download"
        elif basename.startswith("compress_"):
            rule_name = "compress"
        elif basename.startswith("qc_"):
            rule_name = "qc"
        elif basename.startswith("assembly_"):
            rule_name = "assembly"
        elif basename.startswith("rgi_"):
            rule_name = "rgi"
        elif basename.startswith("mobsuite_"):
            rule_name = "mobsuite"
        elif basename.startswith("isescan_"):
            rule_name = "isescan"
        elif basename.startswith("integron_"):
            rule_name = "integron"
        elif basename.startswith("coloc_"):
            rule_name = "coloc"
        elif basename.startswith("network_"):
            rule_name = "network"
        elif basename.startswith("stats_"):
            rule_name = "stats"
        else:
            rule_name = basename.split("_")[0]
            
        try:
            df = pd.read_csv(f, sep="\t")
            if not df.empty:
                row = df.iloc[0]
                records.append({
                    "file": basename,
                    "rule": rule_name,
                    "seconds": float(row["s"]),
                    "max_rss_mb": float(row["max_rss"]),
                    "cpu_time": float(row["cpu_time"]),
                    "mean_load": float(row["mean_load"])
                })
        except Exception as e:
            print(f"Error parsing {f}: {e}")
            
    df_all = pd.DataFrame(records)
    
    # Calculate statistics per rule
    summary = df_all.groupby("rule").agg(
        count=("seconds", "count"),
        mean_seconds=("seconds", "mean"),
        max_seconds=("seconds", "max"),
        mean_rss_mb=("max_rss_mb", "mean"),
        max_rss_mb=("max_rss_mb", "max"),
        mean_cpu_time=("cpu_time", "mean"),
        max_cpu_time=("cpu_time", "max")
    ).reset_index()
    
    # Format times for readability
    summary["mean_time_hms"] = pd.to_timedelta(summary["mean_seconds"], unit="s").astype(str).map(lambda x: x.split(".")[0])
    summary["max_time_hms"] = pd.to_timedelta(summary["max_seconds"], unit="s").astype(str).map(lambda x: x.split(".")[0])
    
    print("=== BENCHMARK SUMMARY PER RULE ===")
    print(summary[["rule", "count", "mean_time_hms", "max_time_hms", "mean_rss_mb", "max_rss_mb"]].to_string(index=False))
    
    # Write to csv
    summary.to_csv("benchmark_summary.csv", index=False)
    print("\nSummary saved to benchmark_summary.csv")

if __name__ == "__main__":
    main()
