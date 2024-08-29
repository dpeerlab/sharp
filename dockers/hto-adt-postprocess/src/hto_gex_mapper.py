import os
import csv
import json
import gzip
import pickle
from tqdm import tqdm

default_path_mapper = "./data/10x-hto-gex-mapper.pickle"


def create(path_10x_whitelist: str) -> dict:

    # create a mapper (TotalSeq-B/C HTO <--> GEX mapper)
    mapper = dict()

    if path_10x_whitelist.endswith(".txt"):
        with open(path_10x_whitelist, "r") as csv_file:
            csv_reader = csv.reader(csv_file, delimiter="\t")
            for row in tqdm(csv_reader, disable=None):
                mapper[row[0].strip()] = row[0].strip()
    elif path_10x_whitelist.endswith(".txt.gz"):
        with gzip.open(path_10x_whitelist, "rt") as csv_file:
            csv_reader = csv.reader(csv_file, delimiter="\t")
            for row in tqdm(csv_reader, disable=None):
                mapper[row[0].strip()] = row[1].strip()
    else:
        raise ValueError(f"File extension not supported: {path_10x_whitelist}")

    return mapper


def write(mapper: dict, path_mapper: str):

    with open(path_mapper, "wb") as fout:
        pickle.dump(mapper, fout)


def load(path_mapper: str) -> dict:

    with open(path_mapper, "rb") as fin:
        mapper = pickle.load(fin)

    return mapper


def exists(path_mapper: str = default_path_mapper) -> bool:

    return os.path.exists(path_mapper)


def decide_which_whitelist(chemistry):
    """
    Based on 10x information, decide which whitelist to use. (https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist-)
    Chemistries are defined in the emulsion sheet.
    Not all whitelists are supported yet.
    """
    with open("./data/whitelists.json") as fin:
        whitelists = json.load(fin)

    for file, chemistries in whitelists.items():
        if chemistry in chemistries:
            return os.path.join("./data", file)
    else:
        raise ValueError("Chemistry {} not supported yet.".format(chemistry))

if __name__ == "__main__":

    mapper_3p = create("./data/3M-february-2018.txt.gz")
    mapper_5p = create("./data/737K-august-2016.txt")

    write(mapper_3p, "./data/10x-hto-gex-mapper-3-prime.pickle")
    write(mapper_5p, "./data/10x-hto-gex-mapper-5-prime.pickle")

    print("DONE.")
