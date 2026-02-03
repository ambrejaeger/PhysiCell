"""Utility functions used by the sensibility analysis scripts.

This module provides helpers to parse PhysiCell XML labels/metadata,
modify XML settings, extract and plot data from output files, combine and
process result files, and run basic Sobol sensitivity analysis tasks.
"""

import os
import glob
import scipy.io
import numpy as np
import xml.etree.ElementTree as ET
import re
import numpy as np
import matplotlib.pyplot as plt
from SALib.analyze.sobol import analyze
from SALib.sample.sobol import sample
import pandas as pd
import shutil
from typing import List, Tuple



def parse_physicell_labels(xml_file:str) -> dict:
    """Parse PhysiCell saved XML labels section and return index ranges and units 
        of the saved parameters in the .mat files.

    :param xml_file: Path to the PhysiCell XML file containing a <labels>
        section.
    :returns: A dictionary mapping label names to a two-element list where
        the first element is a NumPy index array for the label fields and
        the second element is the units string. Returns an empty dict on
        failure.
    """
    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()

        # Find the labels section in the XML
        labels_elem = root.find(".//labels")
        if labels_elem is None:
            raise ValueError(f"Could not find labels in XML file: {xml_file}")

        # Parse each label
        labels = {}
        count = 0
        for label in labels_elem.findall("label"):
            index = int(label.get("index"))
            size = int(label.get("size"))
            units = label.get("units", "none")
            name = label.text.strip() if label.text else f"field_{index}"

            labels[name] = [np.arange(count, count + size), units]
            count += size
        return labels

    except Exception as e:
        print(f"Error parsing XML file {xml_file}: {str(e)}")
        return {}


def parse_physicell_metadata(xml_file:str) -> dict:
    """Extract basic metadata fields from a PhysiCell XML file.

    :param xml_file: Path to the PhysiCell XML file.
    :returns: A dict with entries for `current_time` and `current_runtime`
        if present, otherwise an empty dict on failure.
    """
    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()

        # Find the labels section in the XML
        metadata = root.find(".//metadata")
        if metadata is None:
            raise ValueError(f"Could not find metadata in XML file: {xml_file}")

        # Parse each label
        metadata_dict = {}
        metadata_dict["current_time"] = metadata.find("current_time")
        metadata_dict["current_runtime"] = metadata.find("current_runtime")

        return metadata_dict

    except Exception as e:
        print(f"Error parsing XML file {xml_file}: {str(e)}")
        return {}
    
def update_config_file(config_file, params,values):
    """Update multiple parameters in a PhysiCell XML configuration file.

    :param config_file: Path to the XML configuration file to modify.
    :param params: List of parameter paths to modify.
    :param values: List of new values corresponding to each parameter.
    :returns: True if all modifications were successful, False otherwise.
    """
    all_successful = True
    if len(params) != len(values):
        print("Error: params and values lists must have the same length.")
        return False
    
    for param, value in zip(params, values):
        success = modify_xml(
            xml_file=config_file,
            path=param[0],
            value=value,
            name_cell_def=param[1] if len(param) > 1 else "",
            name_interact_cell_def=param[2] if len(param) > 2 else "",
            variable_name=param[3] if len(param) > 3 else "",
            substrate=param[4] if len(param) > 4 else "",
        )
        if not success:
            print(f"Failed to modify parameter: {param} to value: {value}")
            all_successful = False
    return all_successful

def modify_xml(
    xml_file:str,
    path:str,
    value:str,
    name_cell_def:str="",
    name_interact_cell_def:str="",
    variable_name:str="",
    substrate:str="",
) -> bool:
    """Modify a value in a PhysiCell XML settings file.

    The function supports direct path modification, or targeted updates
    when the path includes markers such as ``cell_definition/`` or
    ``variable/`` and a name is provided.

    :param xml_file: Path to the XML file to modify.
    :param path: XPath-like path to the attribute or element to change.
    :param value: New value to set (will be stringified).

    :param name_cell_def: Optional name of the cell definition to search
        under when using ``cell_definition/`` in the path.
    :param name_interact_cell_def: Optional interaction cell definition name.
    :param variable_name: Optional variable name when using ``variable/`` in
        the path.
    :param substrate: Optional substrate name when searching within a
        cell_definition.
    :returns: True if the file was modified and written, False otherwise.
    """
    try:
        if not os.path.exists(xml_file):
            raise FileNotFoundError(f"XML file '{xml_file}' does not exist")

        target_interaction = None
        subroot = None
        path_to_attribute = ""
        element = None

        modified = False
        tree = ET.parse(xml_file)
        root = tree.getroot()

        if root.find(path) is not None:
            element = root.find(path)

        if variable_name != "":
            path_to_variable = path.split("variable/", 2)[0]
            path_to_attribute = "./" + path.split("variable/", 2)[1]

            if root.find(path_to_variable) is None:
                raise ValueError(
                    f"Path: '{path_to_variable}' not found in XML structure"
                )

            for var in root.findall(os.path.join(path_to_variable, "variable")):
                if var.get("name") == variable_name:
                    subroot = var
                    break

            if subroot is None:
                raise ValueError(
                    f"Variable: '{variable_name}' not found in XML structure"
                )

        if name_cell_def != "":
            path_to_cell_def = path.split("cell_definition/", 2)[0]
            path_to_attribute = "./" + path.split("cell_definition/", 2)[1]

            if root.find(path_to_cell_def) is None:
                raise ValueError(
                    f"Path: '{path_to_cell_def}' not found in XML structure"
                )

            for cell_def in root.findall(
                os.path.join(path_to_cell_def, "cell_definition")
            ):
                if cell_def.get("name") == name_cell_def:
                    subroot = cell_def
                    break
            if subroot is None:
                raise ValueError(
                    f"Cell definition: '{name_cell_def}' not found in XML structure"
                )

            if substrate != "":
                path_to_substrate = "./" + path_to_attribute.split("substrate/", 2)[0]
                path_to_attribute = "./" + path.split("substrate/", 2)[1]

                if subroot.find(path_to_substrate) is None:
                    raise ValueError(
                        f"Path: '{path_to_substrate}' not found in XML structure"
                    )

                for sub in subroot.findall(
                    os.path.join(path_to_substrate, "substrate")
                ):
                    if sub.get("name") == substrate:
                        subroot = sub
                        break

            if name_interact_cell_def != "":
                for interact_cell_def in subroot.findall(path_to_attribute):
                    if interact_cell_def.get("name") == name_interact_cell_def:
                        target_interaction = interact_cell_def
                        break

                if target_interaction is None:
                    error_msg = f"Interaction cell definition '{name_interact_cell_def}' not found in cell definition '{name_cell_def}' at path '{path_to_attribute}'"
                    raise ValueError(error_msg)

        # if (name_cell_def == "") & (name_interact_cell_def == "") & (variable_name == "") & (substrate == ""):

        if target_interaction:
            target_interaction.text = str(value)
            modified = True
        elif subroot != None and (path_to_attribute != ""):
            element = subroot.find(path_to_attribute)

        if element != None:
            element.text = str(value)
            modified = True
        else:
            raise ValueError("Invalid path in the xml")

        if modified:
            tree.write(xml_file)
            return True
        else:
            return False

    except FileNotFoundError as e:
        print(f"Error: {e}")
        return False
    except ET.ParseError as e:
        print(f"Error: Invalid XML format in file '{xml_file}': {e}")
        return False
    except ValueError as e:
        print(f"Error: {e}")
        return False
    except Exception as e:
        print(f"Unexpected error: {e}")
        return False

def get_cells_data(labels_file: str, mat_file: str) -> pd.DataFrame:
    """Load cells data from a .mat file and return a DataFrame indexed by feature.

    In the PhysiCell MAT output the rows correspond to features (position,
    cell_type, etc.) and columns to individual cells.

    :param labels_file: Path to the PhysiCell XML labels file.
    :param mat_file: Path to the MATLAB .mat file containing the `cells`
        array.
    :returns: A pandas DataFrame with feature rows (index) and one column
        per cell.
    """
    labels = parse_physicell_labels(labels_file)
    if not labels:
        raise ValueError(f"Could not parse labels from: {labels_file}")

    mat = scipy.io.loadmat(mat_file)
    cells_data = mat["cells"]  # expected shape: (n_features, n_cells)

    from collections import defaultdict

    name_counts = defaultdict(int)
    feature_names = []
    row_arrays = []

    for label_name, (indices, _units) in labels.items():
        idxs = np.atleast_1d(indices).astype(int).ravel()
        for idx in idxs:
            name_counts[label_name] += 1
            suffix = f"_{name_counts[label_name]}" if len(idxs) > 1 else ""
            feat_name = f"{label_name}{suffix}"

            try:
                row = np.asarray(cells_data[idx, :]).ravel()
            except Exception as e:
                raise IndexError(f"Error reading feature index {idx} from '{mat_file}': {e}")

            feature_names.append(feat_name)
            row_arrays.append(row)

    if not row_arrays:
        return pd.DataFrame()

    data_matrix = np.vstack(row_arrays)  # shape (n_features, n_cells)
    n_cells = data_matrix.shape[1]

    idx = labels['ID'][0].astype(int)
    col_names = [f"cell_{int(cells_data[idx, i].item())}" for i in range(n_cells)]

    df = pd.DataFrame(data_matrix, index=feature_names, columns=col_names)
    return df

def extract_position_type_data(mat_file: str, xml_file: str) -> tuple:
    """Load a .mat cells file and extract positions, types and IDs.

    :param mat_file: Path to the MATLAB .mat file produced by PhysiCell.
    :param xml_file: Path to the corresponding PhysiCell XML file used to
        interpret the columns (labels section).
    :type xml_file: str
    :returns: Tuple ``(positions, cells_type, cells_ID)`` where ``positions``
        is an (N, 3) NumPy array and ``cells_type``/``cells_ID`` are 1-D
        arrays. If an error occurs the function returns ``(None, None, None)``.
    :rtype: tuple
    """
    try:
        labels = parse_physicell_labels(xml_file)
        if not labels:
            error_msg = (
                f"Impossible to parse_physicell_labels of file at path {xml_file}"
            )
            raise ValueError(error_msg)

        # Load mat file
        data = scipy.io.loadmat(mat_file)
        cells_data = data["cells"]

        print(f"Processing {os.path.basename(mat_file)}")
        print(f"Cells array shape: {cells_data.shape}")

        num_cells = cells_data.shape[1]

        positions = np.empty((num_cells, 3))
        cells_type = np.empty(num_cells)
        cells_ID = np.empty(num_cells)

        position_indices = labels["position"][0]
        id_index = labels["ID"][0]
        cell_type_index = labels["cell_type"][0]

        for i in range(num_cells):
            try:
                positions[i] = cells_data[position_indices, i]
                cells_type[i] = cells_data[cell_type_index, i].item()
                cells_ID[i] = cells_data[id_index, i].item()

            except (IndexError, ValueError) as e:
                print(f"Warning: Error processing cell {i}: {e}")
                continue

        return positions, cells_type, cells_ID

    except Exception as e:
        print(f"Error processing file {os.path.basename(mat_file)}: {str(e)}")
        return None, None, None


def plot_cells_2D(positions):
    mask_x0 = positions[:, 0] != 0.0
    x = positions[mask_x0, 0]

    y = positions[mask_x0, 1]
    plt.figure(figsize=(12, 6))
    print("y:", y)
    # Plot line
    plt.plot(x, y, "b-", linewidth=3, label="data")

    # Customize the plot
    plt.xlabel("X values", fontsize=12)
    plt.ylabel("Y values", fontsize=12)
    plt.grid(True, alpha=0.3, linestyle="--")
    plt.legend(fontsize=10)
    plt.tight_layout()

    plt.show()


def get_output_files(path, prefix="output", suffix="_cells.mat"):

    pattern = os.path.join(path, prefix + "*" + suffix)
    files = glob.glob(pattern)

    pattern_regex = re.compile(rf"{re.escape(prefix)}(\d+){re.escape(suffix)}$")

    def extract_number(filepath):
        basename = os.path.basename(filepath)
        match = pattern_regex.search(basename)
        if match:
            return int(match.group(1))
        return -1

    files.sort(key=extract_number)

    return files


def define_problem(num_vars, names, bounds, groups=[]):
    if len(groups) == 0:
        problem = {"num_vars": num_vars, "names": names, "bounds": bounds}
    else:
        problem = {
            "groups": groups,
            "num_vars": num_vars,
            "names": names,
            "bounds": bounds,
        }
    return problem


def define_set_param(num_vars, names, bounds, groups=[], sample_size=32, seed=0):
    print("Defining parameter space and sampling...")
    problem = define_problem(num_vars, names, bounds, groups)

    param_values = sample(
        problem, sample_size, seed=seed
    )  # Génère N*(2+D) jeux de paramètres avec D le nombre de paramètres et N un multiple de 2 fourni en argument
    return param_values


def analyze_sobol(result_file, param_names_file, bounds, groups=[], column=1):
    import ast

    with open(param_names_file, "r") as f:
        names = [ast.literal_eval(line.strip()) for line in f if line.strip()]
    names = np.asarray([f"{';'.join(map(str, sublist))}" for sublist in names])

    problem = define_problem(len(names), names, bounds, groups)
    print(problem)

    with open(result_file, "r") as f:
        lines = f.readlines()

    sorted_lines = sorted(lines, key=lambda x: int(x.split()[0]))
    Y = np.asarray([float(line.split()[column]) for line in sorted_lines])
    print("Array:", Y)
    print("Length:", len(Y))

    Si = analyze(problem, Y, print_to_console=True)
    return Si


###########  Functions for output processing and analysis ###########

def process_file(input_file: str, output_file: str) -> List[str]:
    """
    Suppresses lines with duplicate first columns, orders the output by the first column,
    and returns the deleted lines.
    
    :param input_file: Path to the input text file.
    :param output_file: Path where the sorted, unique-key file will be written.
    :returns: A list of lines that were suppressed/deleted.
    """
    seen_keys = set()
    kept_entries = []
    deleted_lines = []

    with open(input_file, "r") as f:
        for line in f:
            clean_line = line.strip()
            if not clean_line:
                continue

            parts = clean_line.split()
            if not parts:
                continue
            
            # The 'key' is the first column
            first_col_raw = parts[0]
            
            if first_col_raw in seen_keys:
                deleted_lines.append(clean_line)
            else:
                seen_keys.add(first_col_raw)
                # We store the raw string and a sort-friendly key
                try:
                    sort_key = float(first_col_raw)
                except ValueError:
                    sort_key = first_col_raw
                
                kept_entries.append((sort_key, clean_line))

    # Sort the kept entries by the first column (increasing order)
    kept_entries.sort(key=lambda x: x[0])

    # Write the sorted, unique entries to the output file
    with open(output_file, "w") as f:
        for _, line_text in kept_entries:
            f.write(line_text + "\n")

    # Reporting
    print(f"--- Process Complete ---")
    print(f"Lines suppressed: {len(deleted_lines)}")
    print(f"Unique lines saved: {len(kept_entries)}")
    print(f"Output saved to: {output_file}\n")

    return deleted_lines

def combine_files_with_header(
    file1_path:str, file2_path:str, output_path:str, column_names:list, column_1:str="nbr_breaks"
) -> pd.DataFrame:
    """Combine two files by assigning values from file1 into rows of file2.

    file1 is expected to contain row indices and values which will be placed
    into a new leading column in the combined output. File2 contains the
    rest of the tabular data.

    :param file1_path: Path to the first file (row indices + values).
    :param file2_path: Path to the second file (tabular data).
    :param output_path: Output path for the combined file.
    :param column_names: Column names to apply to file2 if provided.
    :param column_1: Name for the inserted column from file1 (default
        'nbr_breaks').
    :returns: The combined DataFrame.
    """

    df1 = pd.read_csv(file1_path, sep=r"\s+", header=None, names=["row_idx", column_1])
    df2 = pd.read_csv(file2_path, sep=r"\s+", header=None)

    if len(column_names) == df2.shape[1]:
        df2.columns = column_names

    else:
        print(
            f"Warning: column_names has {len(column_names)} items, but file2 has {df2.shape[1]} columns"
        )
        df2.columns = [f"col_{i+1}" for i in range(df2.shape[1])]

    nbr_breaks_list = [0.0] * len(df2)

    for idx, row in df1.iterrows():
        row_idx = int(row["row_idx"])
        if 0 <= row_idx < len(df2):
            nbr_breaks_list[row_idx] = row[column_1]

    df_combined = pd.DataFrame({column_1: nbr_breaks_list})
    df_combined = pd.concat([df_combined, df2.reset_index(drop=True)], axis=1)

    df_combined.to_csv(output_path, sep="\t", index=False, float_format="%.6e")

    print(f"Combined file created: {output_path}")
    print(f"Shape of combined data: {df_combined.shape}")
    print(f"Rows from file1 assigned: {len(df1[df1['row_idx'] < len(df2)])}")

    return df_combined


def combine_files_with_header_2(
    file1_path:str, file2_path:str, output_path:str, column_names:list, column_1:str, 
    column_2:str) -> pd.DataFrame:

    """Combine two files where file1 provides two named columns to insert.

    Similar to :func:`combine_files_with_header` but extracts two value
    columns from file1 (specified by column_1 and column_2) and inserts
    them into the combined output.

    :param file1_path: Path to the first file.
    :param file2_path: Path to the second file.
    :param output_path: Path where the combined file will be saved.
    :param column_names: Column names for file2.
    :param column_1: Column name in file1 for the first inserted value.
    :param column_2: Column name in file1 for the second inserted value.
    :returns: The combined DataFrame.
    """

    df1 = pd.read_csv(
        file1_path, sep=r"\s+", header=None, names=["row_idx", column_1, column_2]
    )
    df2 = pd.read_csv(file2_path, sep=r"\s+", header=None)

    print(df1)

    if len(column_names) == df2.shape[1]:
        df2.columns = column_names

    else:
        print(
            f"Warning: column_names has {len(column_names)} items, but file2 has {df2.shape[1]} columns"
        )
        df2.columns = [f"col_{i+1}" for i in range(df2.shape[1])]

    value_column1_list = [0.0] * len(df2)
    value_common2_list = [0.0] * len(df2)

    for idx, row in df1.iterrows():
        row_idx = int(row["row_idx"])
        if 0 <= row_idx < len(df2):
            value_column1_list[row_idx] = row[column_1]
            value_common2_list[row_idx] = row[column_2]

    df_combined = pd.DataFrame(
        {column_1: value_column1_list, column_2: value_common2_list}
    )
    df_combined = pd.concat([df_combined, df2.reset_index(drop=True)], axis=1)

    df_combined.to_csv(output_path, sep="\t", index=False, float_format="%.6e")

    print(f"Combined file created: {output_path}")
    print(f"Shape of combined data: {df_combined.shape}")

    return df_combined


def extract_X_Y(file_path, X_columns, Y_column="nbr_breaks"):
    """Extract feature matrix X and target Y from a tab-separated file.

    X_columns may be a single column name/index or a list of names/indices.

    :param file_path: Path to the tab-separated file.
    :type file_path: str
    :param X_columns: Column(s) to use as features (str, int, or list).
    :type X_columns: str|int|list
    :param Y_column: Column name for the target variable (default
        'nbr_breaks').
    :type Y_column: str
    :returns: Tuple (X, Y) where X is a list or list-of-lists and Y is a 1-D
        NumPy array.
    :rtype: tuple
    """

    df = pd.read_csv(file_path, sep="\t")
    if Y_column not in df.columns:
        raise ValueError(
            f"Y column '{Y_column}' not found in file. Available columns: {list(df.columns)}"
        )

    Y = df[Y_column].values
    if isinstance(X_columns, str):
        if X_columns not in df.columns:
            raise ValueError(f"Column '{X_columns}' not found.")
        X = df[X_columns].values.tolist()

    elif isinstance(X_columns, int):
        if X_columns >= len(df.columns) or X_columns < 0:
            raise ValueError(f"Column index {X_columns} out of range.")
        X = df.iloc[:, X_columns].values.tolist()

    elif isinstance(X_columns, list):
        if not X_columns:
            raise ValueError("X_columns list is empty.")

        if all(isinstance(x, str) for x in X_columns):
            missing_cols = [col for col in X_columns if col not in df.columns]
            if missing_cols:
                raise ValueError(f"Columns not found: {missing_cols}")

            if len(X_columns) == 1:
                X = df[X_columns[0]].values.tolist()
            else:
                X = df[X_columns].values.T.tolist()

        elif all(isinstance(x, int) for x in X_columns):
            valid_indices = [idx for idx in X_columns if 0 <= idx < len(df.columns)]
            if len(valid_indices) != len(X_columns):
                raise ValueError(
                    f"Some indices are out of range. Valid range: 0-{len(df.columns)-1}"
                )

            if len(X_columns) == 1:
                X = df.iloc[:, X_columns[0]].values.tolist()
            else:
                X = df.iloc[:, X_columns].values.T.tolist()

        else:
            raise ValueError("X_columns list must contain all strings or all integers.")

    else:
        raise TypeError(
            "X_columns must be string, integer, or list of strings/integers."
        )

    print(f"Successfully extracted:")
    print(f"  Y shape: {len(Y)} values")
    print(f"  X shape: {len(X) if isinstance(X[0], list) else '1D'} features")

    return X, Y


def plot_scatter_sets(
    Y,
    X,
    set_names=None,
    title="Scatter Plot",
    xlabel="Features",
    ylabel="nbr_breaks",
    figsize=(10, 6),
    save_path=None,
    show_plot=True,
):

    """Scatter-plot one or more feature sets against the target Y.

    :param Y: Target values (1-D iterable).
    :type Y: array-like
    :param X: A list (or list of lists) of feature values to plot against Y.
    :type X: list
    :param set_names: Optional list of names for each feature set.
    :type set_names: list
    :param title: Plot title.
    :type title: str
    :param xlabel: X-axis label.
    :type xlabel: str
    :param ylabel: Y-axis label.
    :type ylabel: str
    :param figsize: Figure size tuple.
    :type figsize: tuple
    :param save_path: Optional path to save the figure.
    :type save_path: str
    :param show_plot: If True, display the plot interactively.
    :type show_plot: bool
    :returns: Matplotlib Figure and Axes objects.
    :rtype: tuple
    """

    if isinstance(X, list):
        if X and isinstance(X[0], list):
            num_sets = len(X)
            for i, x_set in enumerate(X):
                if len(x_set) != len(Y):
                    raise ValueError(
                        f"Set {i} has {len(x_set)} values, but Y has {len(Y)} values"
                    )
        else:
            num_sets = 1
            X = [X]
    else:
        raise TypeError("X must be a list or list of lists")

    fig, ax = plt.subplots(figsize=figsize)

    # Default set names if not provided
    if set_names is None:
        set_names = [f"Feature Set {i+1}" for i in range(num_sets)]
    elif len(set_names) != num_sets:
        print(
            f"Warning: {len(set_names)} names provided for {num_sets} sets. Using default names."
        )
        set_names = [f"Feature Set {i+1}" for i in range(num_sets)]

    for i in range(num_sets):
        x_values = X[i]
        ax.scatter(
            x_values, Y, alpha=0.6, edgecolors="w", linewidth=0.5, label=set_names[i]
        )

    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Plot saved to: {save_path}")

    if show_plot:
        plt.show()

    return fig, ax


def save_dataframes_to_txt(df_list, filename, df_names=None, include_index=True):
    """Save a list of pandas DataFrames to a human-readable text file.

    Each DataFrame is written as a section beginning with ``=== DataFrame: ...``.

    :param df_list: List of pandas DataFrame objects to save.
    :type df_list: list
    :param filename: Path to the output text file.
    :type filename: str
    :param df_names: Optional list of names for each DataFrame.
    :type df_names: list
    :param include_index: Whether to include the DataFrame index in output.
    :type include_index: bool
    :returns: None
    """
    with open(filename, "w") as f:
        for i, df in enumerate(df_list):
            if df_names and i < len(df_names):
                f.write(f"=== DataFrame: {df_names[i]} ===\n")
            else:
                f.write(f"=== DataFrame {i+1} ===\n")

            df_str = df.to_string(index=include_index, header=True)
            f.write(df_str)
            f.write("\n\n")

    print(f"Saved {len(df_list)} DataFrames to '{filename}'")


def load_dataframes_from_txt(filename):
    """Load DataFrames previously saved with :func:`save_dataframes_to_txt`.

    The loader attempts to parse each ``=== DataFrame ...`` section back
    into a pandas DataFrame. Sections that cannot be parsed are skipped
    with a warning.

    :param filename: Path to the input text file.
    :type filename: str
    :returns: List of reconstructed pandas DataFrame objects.
    :rtype: list
    """
    dataframes = []
    current_df_lines = []
    reading_df = False

    with open(filename, "r") as f:
        lines = f.readlines()

    for line in lines:
        # Check if DataFrame header
        if line.startswith("=== DataFrame"):
            if reading_df and current_df_lines:
                df_str = "".join(current_df_lines)
                try:
                    lines_split = df_str.strip().split("\n")
                    if len(lines_split) > 1:
                        df = pd.read_csv(
                            pd.io.common.StringIO(df_str), sep=r"\s+", engine="python"
                        )
                        dataframes.append(df)
                except:
                    print(f"Warning: Could not parse DataFrame section")

                current_df_lines = []

            reading_df = True
            continue

        if reading_df and line.strip():
            current_df_lines.append(line)

    if reading_df and current_df_lines:
        df_str = "".join(current_df_lines)
        try:
            df = pd.read_csv(pd.io.common.StringIO(df_str), sep=r"\s+", engine="python")
            dataframes.append(df)
        except:
            print(f"Warning: Could not parse last DataFrame section")

    return dataframes

if __name__ == "__main__":
    labels_file = "./output/initial.xml"
    mat_file = "./output/first_save_cells.mat"
    data = get_cells_data(labels_file, mat_file)
    print(data)
