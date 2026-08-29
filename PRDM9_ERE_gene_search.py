# Search for PRDM9 & ERE genes in HOMER output
## CONCLUSION: no significant results

import os
import re
import sys

BASE = r"PATH TO HOMER DATA FOLDER"   # replace with path

WAGE_KNOWN   = os.path.join(BASE, "DMRs_wAGE", "motifs_wAGE_output", "knownResults_wAGE")
WAGE_DENOVO  = os.path.join(BASE, "DMRs_wAGE", "motifs_wAGE_output", "homerResults_wAGE")
NOAGE_KNOWN  = os.path.join(BASE, "DMRs_woAGE", "motif_woAGE_output2", "knownResults")
NOAGE_DENOVO = os.path.join(BASE, "DMRs_woAGE", "motif_woAGE_output2", "homerResults")

WAGE_KNOWN_TXT    = os.path.join(BASE, "DMRs_wAGE", "motifs_wAGE_output", "knownResults_wAGE.txt")
WAGE_KNOWN_HTML   = os.path.join(BASE, "DMRs_wAGE", "motifs_wAGE_output", "knownResults_wAGE.html")
WAGE_DENOVO_HTML  = os.path.join(BASE, "DMRs_wAGE", "motifs_wAGE_output", "homerResults_wAGE.html")
NOAGE_KNOWN_TXT   = os.path.join(BASE, "DMRs_woAGE", "motif_woAGE_output2", "knownResults.txt")
NOAGE_KNOWN_HTML  = os.path.join(BASE, "DMRs_woAGE", "motif_woAGE_output2", "knownResults.html")
NOAGE_DENOVO_HTML = os.path.join(BASE, "DMRs_woAGE", "motif_woAGE_output2", "homerResults.html")

# =========================================
# Search terms
# =========================================
PRDM9_TERMS = ["PRDM9"]
ERE_TERMS   = ["estrogen", "ERE", "ESR1", "ESR2", "ER-alpha", "ER-beta", "NR3C"]

# =========================================
# Output file
# =========================================
OUTPUT_FILE = os.path.join(BASE, "motif_search_results.txt")

# =========================================
# Helper functions
# =========================================
def search_file(filepath, terms):
    if not os.path.exists(filepath):
        return None
    matches = []
    pattern = re.compile("|".join(terms), re.IGNORECASE)
    with open(filepath, "r", errors="ignore") as f:
        for line_num, line in enumerate(f, 1):
            if pattern.search(line):
                matches.append((line_num, line.strip()))
    return matches

def search_and_report(label, filepaths, terms):
    print(f"\n--- {label} ---")
    any_found = False
    for filepath in filepaths:
        filename = os.path.basename(filepath)
        results = search_file(filepath, terms)
        if results is None:
            print(f"  [{filename}]: File not found")
        elif len(results) == 0:
            print(f"  [{filename}]: No matches found")
        else:
            any_found = True
            print(f"  [{filename}]: {len(results)} match(es) found")
            for line_num, line in results:
                display = line[:300] + "..." if len(line) > 300 else line
                print(f"    Line {line_num}: {display}")
    return any_found

# =========================================
# Redirect all print output to file AND terminal
# =========================================
class Tee:
    """Writes output to both terminal and a file simultaneously."""
    def __init__(self, filepath):
        self.terminal = sys.stdout
        self.file = open(filepath, "w", encoding="utf-8")

    def write(self, message):
        self.terminal.write(message)
        self.file.write(message)

    def flush(self):
        self.terminal.flush()
        self.file.flush()

    def close(self):
        self.file.close()

# Start capturing output
tee = Tee(OUTPUT_FILE)
sys.stdout = tee

# =========================================
# Run searches
# =========================================
print("=========================================")
print("Searching HOMER Results for PRDM9 and ERE")
print("=========================================")

print("\n\n=========================================")
print("PRDM9 RESULTS")
print("=========================================")
search_and_report("With AGE: Known Results",    [WAGE_KNOWN_TXT, WAGE_KNOWN_HTML],   PRDM9_TERMS)
search_and_report("With AGE: De Novo Results",  [WAGE_DENOVO_HTML],                  PRDM9_TERMS)
search_and_report("Without AGE: Known Results", [NOAGE_KNOWN_TXT, NOAGE_KNOWN_HTML], PRDM9_TERMS)
search_and_report("Without AGE: De Novo",       [NOAGE_DENOVO_HTML],                 PRDM9_TERMS)

print("\n\n=========================================")
print("ESTROGEN RESPONSE ELEMENT RESULTS")
print("=========================================")
search_and_report("With AGE: Known Results",    [WAGE_KNOWN_TXT, WAGE_KNOWN_HTML],   ERE_TERMS)
search_and_report("With AGE: De Novo Results",  [WAGE_DENOVO_HTML],                  ERE_TERMS)
search_and_report("Without AGE: Known Results", [NOAGE_KNOWN_TXT, NOAGE_KNOWN_HTML], ERE_TERMS)
search_and_report("Without AGE: De Novo",       [NOAGE_DENOVO_HTML],                 ERE_TERMS)

print("\n\n=========================================")
print("Search Complete")
print("=========================================")

# =========================================
# Close file and restore terminal output
# =========================================
sys.stdout = tee.terminal
tee.close()

print(f"\nResults saved to: {OUTPUT_FILE}")