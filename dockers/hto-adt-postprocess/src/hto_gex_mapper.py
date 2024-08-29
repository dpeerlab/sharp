import os
import json

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

