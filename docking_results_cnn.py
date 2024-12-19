import os
import pandas as pd
import re

folder_path = "/content/drive/MyDrive/12-13.12 results"

# Function to parse a .txt file and select the best result
def analyze_txt_file(file_path):
    best_pose = None
    results = []

    with open(file_path, "r") as file:
        lines = file.readlines()

        # Searching for the beginning of the table
        is_table = False
        for line in lines:
            if line.strip().startswith("-----+"):
                is_table = True
                continue

            if is_table and re.match(r"\s*\d+", line):  # Lines starting with a digit
                try:
                    # Extracting data from the line using a regular expression
                    match = re.match(
                        r"\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+.+\s+(-?\d+\.\d+)", line
                    )
                    if match:
                        mode = int(match.group(1))        # First column - position number
                        affinity = float(match.group(2))  # Second column - affinity score
                        intramol = float(match.group(3))  # Third column - intramolecular
                        cnn_affinity = float(match.group(4))  # Fifth column - CNN affinity
                        results.append({
                            "mode": mode,
                            "affinity": affinity,
                            "intramol": intramol,
                            "cnn_affinity": cnn_affinity
                        })
                except (ValueError, IndexError):
                    continue

    # Selecting the best result (lowest affinity, highest intramolecular, lowest CNN affinity)
    if results:
        best_pose = min(results, key=lambda x: (x["cnn_affinity"], x["affinity"], -x["intramol"]))

    return best_pose

# Analyzing all .txt files in the folder
best_results = []
for file_name in os.listdir(folder_path):
    if file_name.endswith(".txt"):
        file_path = os.path.join(folder_path, file_name)
        best_pose = analyze_txt_file(file_path)
        if best_pose:
            best_pose["file"] = file_name
            best_results.append(best_pose)

# Creating the global ranking
if best_results:
    global_ranking = pd.DataFrame(best_results)
    global_ranking = global_ranking.sort_values(
        by=["cnn_affinity", "affinity", "intramol"],
        ascending=[True, True, False] 
    )

    # Saving the ranking to a CSV file
    global_ranking.to_csv("global_ranking_cnn.csv", index=False)

    print("The global ranking has been saved to 'global_ranking_cnn.csv'.")
else:
    print("No results were found in any .txt file.")
