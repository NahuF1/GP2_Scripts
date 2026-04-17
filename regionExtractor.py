import os
import pathlib
import csv
import subprocess
from datetime import datetime
from multiprocessing import Manager
from concurrent.futures import ProcessPoolExecutor, as_completed

release = 11
RELEASE_PATH = pathlib.Path(pathlib.Path.home(), f'workspace/gp2_tier2_eu_release{release}')
dataset = ["NBA", "WGS"]
PATH_NBA_GENO = pathlib.Path(RELEASE_PATH, 'imputed_genotypes')
PATH_WGS_GENO = pathlib.Path(RELEASE_PATH, 'wgs/deepvariant_joint_calling/plink')
ANCESTRIES = ['AAC', 'AFR', 'AJ', 'AMR', 'CAS', 'EAS', 'EUR', 'FIN', 'MDE', 'SAS', 'CAH']
#
DIR_TOOL = "/home/jupyter/tools"
DIR_WSPS = "/home/jupyter/workspace/ws_files"
DIR_ORIG = "/home/jupyter/workspace/ws_files/Original_Files"

with open(f'{DIR_ORIG}/HAR_list_phase_1.tsv', 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    HARS_DICT = {
        row[3]: {
            'name':  row[3],
            'chrom': row[0].replace('chr', ''),
            'start': int(row[1]),
            'end':   int(row[2]),
        }
        for row in reader
    }

def regionExtractor(HAR, chrom, startBP, endBP, sets, PATH_GENO, ANCESTRY, exceptions_file, lock):
    MAIN = f"/home/jupyter/workspace/ws_files/Working_{sets}/{ANCESTRY}/InputFiles"

    ts_start = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[START] {HAR} - {ANCESTRY} | chr{chrom}:{startBP}:{endBP} | {ts_start}")

    inputRawFile = f"{PATH_GENO}/{ANCESTRY}/chr{chrom}_{ANCESTRY}_release11_vwb"
    outputDir = f"{MAIN}/Indiv_HARS"
    covar = f"{MAIN}/{ANCESTRY}_covariate_file.txt"
    outputPrefix = f"{outputDir}/{HAR}"
    os.makedirs(outputDir, exist_ok=True)

    regionExtractorNBA = [
        f"{DIR_TOOL}/plink2",
        "--pfile", inputRawFile,
        "--chr", str(chrom),
        "--from-bp", str(startBP),
        "--to-bp", str(endBP),
        "--extract-if-info", "R2>=0.8",
        "--mac", "2",
        "--hwe", "0.0001", "keep-fewhet",
        "--max-maf", "0.05",
        "--make-pgen",
        "--out", outputPrefix
    ]
    regionExtractorWGS = [
        f"{DIR_TOOL}/plink2",
        "--pfile", inputRawFile,
        "--chr", str(chrom),
        "--from-bp", str(startBP),
        "--to-bp", str(endBP),
        "--mac", "2",
        "--max-maf", "0.05",
        "--make-pgen",
        "--pheno", covar,
        "--not-pheno", "FATID", "MATID", "SEX", "AGE", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10",
        "--out", outputPrefix
    ]

    if sets == "NBA":
        result = subprocess.run(regionExtractorNBA, capture_output=True, text=True)
    elif sets == "WGS":
        result = subprocess.run(regionExtractorWGS, capture_output=True, text=True)

    if "No variants remaining after main filters" in result.stderr:
        with lock:
            with open(exceptions_file, 'a') as f:
                f.write(f'{HAR}\t{ANCESTRY}\t{chrom}\t{startBP}\t{endBP}\n')
        return f"[SKIP]  {HAR} - {ANCESTRY} | no variants passed filters"

    else:
        try:
            vcfConverter = [
                f"{DIR_TOOL}/plink2",
                "--pfile", outputPrefix,
                "--recode", "vcf", "id-paste=iid",
                "--out", outputPrefix
            ]
            subprocess.run(vcfConverter, capture_output=True, check=True)
            subprocess.run(["bgzip", "-f", f"{outputPrefix}.vcf"], check=True, capture_output=True)
            subprocess.run(["tabix", "-f", "-p", "vcf", f"{outputPrefix}.vcf.gz"], check=True, capture_output=True)

        except subprocess.CalledProcessError as e:
            return f"[ERROR] {HAR} - {ANCESTRY} | VCF conversion failed: {e}"

    ts_end = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[END]   {HAR} - {ANCESTRY} | chr{chrom}:{startBP}:{endBP} | {ts_end}")
    return f"[DONE]  {HAR} - {ANCESTRY} | chr{chrom}:{startBP}:{endBP}"


def run_task(args):
    HAR, chrom, startBP, endBP, sets, PATH_GENO, ANCESTRY, exceptions_file, lock = args
    return regionExtractor(HAR, chrom, startBP, endBP, sets, PATH_GENO, ANCESTRY, exceptions_file, lock)


if __name__ == "__main__":
    BUILD = "hg38"
    max_workers = 28

    with Manager() as manager:
        lock = manager.Lock()

        for sets in dataset:
            PATH_GENO = str(PATH_NBA_GENO if sets == "NBA" else PATH_WGS_GENO)

            # Pre-create directories and exception files
            for ANCESTRY in ANCESTRIES:
                exceptions_file = str(pathlib.Path(DIR_WSPS, f'Working_{sets}', ANCESTRY, f'failed_HARs_{sets}_{ANCESTRY}.tsv'))
                with open(exceptions_file, 'w') as f:
                    f.write('HAR\tancestry\tchr\tstart\tend\n')
                outDir = pathlib.Path(DIR_WSPS, f'Working_{sets}', ANCESTRY, 'InputFiles', 'Indiv_HARS')
                os.makedirs(outDir, exist_ok=True)

            # Build task list
            tasks = [
                (HAR,
                 HARS_DICT[HAR]["chrom"],
                 str(HARS_DICT[HAR]["start"]),
                 str(HARS_DICT[HAR]["end"]),
                 sets,
                 PATH_GENO,
                 ANCESTRY,
                 str(pathlib.Path(DIR_WSPS, f'Working_{sets}', ANCESTRY, f'failed_HARs_{sets}_{ANCESTRY}.tsv')),
                 lock)
                for ANCESTRY in ANCESTRIES
                for HAR in HARS_DICT
            ]

            print(f"\nTotal tasks: {len(tasks)} ({len(HARS_DICT)} HARs × {len(ANCESTRIES)} ancestries)")
            print(f"Running with {max_workers} workers\n")

            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                futures = {executor.submit(run_task, task): task for task in tasks}

                completed = 0
                for future in as_completed(futures):
                    task = futures[future]
                    completed += 1
                    try:
                        result = future.result()
                        print(f"[{completed}/{len(tasks)}] {result}")
                    except Exception as e:
                        print(f"[{completed}/{len(tasks)}] Task failed for {task[0]} - {task[6]}: {e}")