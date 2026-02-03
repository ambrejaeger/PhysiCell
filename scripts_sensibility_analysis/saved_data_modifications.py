"""Helpers to modify and extract saved data.

This module provides utilities for modifying saved .txt files generated from using the 
PhysiSandS addon to PhysiCell and extracting values from PhysiCell XML settings files.

"""

import re
import os
import xml.etree.ElementTree as ET


def modify_cell_data(file_path :str, condition_key: str, condition_val: str, target_path: str, new_value: str) -> bool:
    """Modify cell data in a custom text format based on a condition.

    :param file_path: Path to the .txt file.
    :param condition_key: The key to check within each cell block (e.g. 'type').
    :param condition_val: The value the key must match. Compared as string
        against the file contents.
    :param target_path: Slash-separated path to the target field inside a
        cell block, for example ``"Section/SubSection/Key"``.
    :param new_value: The value to write for the target key.
    :returns: True if the file was modified and saved, False otherwise
        (including when the file does not exist).
    :rtype: bool
    """
    if not os.path.exists(file_path):
        return False

    with open(file_path, "r") as f:
        content = f.read()

    # Split the file into individual Cell blocks
    # We use a lookahead to split by 'Cell' but keep the word 'Cell'
    blocks = re.split(r"(^Cell\n)", content, flags=re.MULTILINE)

    # Reassemble blocks (split creates empty entries and separates 'Cell' from its body)
    refined_blocks = []
    for i in range(1, len(blocks), 2):
        refined_blocks.append(blocks[i] + blocks[i + 1])

    modifications_made = False

    for b_idx, block in enumerate(refined_blocks):
        lines = block.splitlines()

        # 1. Check if this cell meets the condition
        # Look for the condition_key: condition_val
        cell_matches = False
        for line in lines:
            if f"{condition_key}:" in line:
                val = line.split(":")[-1].strip()
                if val == str(condition_val):
                    cell_matches = True
                    print("Cell modified")
                    break

        if not cell_matches:
            continue

        # 2. If matched, navigate the target path
        # Target path format: "Section/SubSection/Key"
        path_parts = target_path.split("/")
        target_key = path_parts[-1]
        sections = path_parts[:-1]

        current_line_idx = 0
        search_range_start = 0
        search_range_end = len(lines)

        # Narrow down the search by iterating through sections
        for section in sections:
            found_section = False
            for i in range(search_range_start, len(lines)):
                if section in lines[i]:
                    search_range_start = i
                    found_section = True
                    break
            if not found_section:
                break

        # 3. Find and replace the value within the found range
        for i in range(search_range_start, len(lines)):
            # If we hit another major section or the end of the cell, stop
            if (
                i > search_range_start
                and ":" in lines[i]
                and not any(s in lines[i] for s in [target_key] + sections)
            ):
                # This is a simple heuristic to ensure we don't jump into the next block section
                pass

            if f"{target_key}:" in lines[i]:
                lines[i] = f"{target_key}: {new_value}"
                modifications_made = True
                break

        refined_blocks[b_idx] = "\n".join(lines) + "\n"

    # Save the file if changes were made
    if modifications_made:
        with open(file_path, "w") as f:
            f.write("".join(refined_blocks))
        return True

    return False


def extract_data_xml(
    xml_file:str,
    path:str,
    name_cell_def:str ="",
    name_interact_cell_def:str ="",
    variable_name:str ="",
    substrate:str ="",
) -> str | None:
    """Extract a value from a PhysiCell XML settings file.

    The function supports several lookup modes: direct path lookup, searching
    for a named variable under a ``variable/`` section, or searching for a
    named cell definition under a ``cell_definition/`` section. It returns
    the text content of the located element or ``False`` on error.

    :param xml_file: Path to the XML file.
    :param path: XPath-like path used to locate the element.
    :param name_cell_def: Name of the cell definition to search under.
    :param name_interact_cell_def: Name of the cell type involved in the 
        interaction to find under the selected cell definition.
    :param variable_name: Name of the variable when the path points into
        ``variable/``.
    :param substrate: Substrate name when searching inside a
        ``cell_definition``.
    :returns: The text content of the found XML element on success, or
        ``None`` if an error occurred.
    :rtype: str or None

    .. note::

        The function will catch common errors (missing file, parse errors,
        and lookup failures), print a short descriptive message and return
        ``None``.
    """
    try:
        if not os.path.exists(xml_file):
            raise FileNotFoundError(f"XML file '{xml_file}' does not exist")

        target_interaction = None
        subroot = None
        path_to_attribute = ""
        element = None

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

        if target_interaction:
            return target_interaction.text
        elif subroot != None and (path_to_attribute != ""):
            element = subroot.find(path_to_attribute)

        if element != None:
            return element.text
        else:
            raise ValueError("Invalid path in the xml")
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return None
    except ET.ParseError as e:
        print(f"Error: Invalid XML format in file '{xml_file}': {e}")
        return None
    except ValueError as e:
        print(f"Error: {e}")
        return None
    except Exception as e:
        print(f"Unexpected error: {e}")
        return None


if __name__ == "__main__":
    print(
        extract_data_xml(
            "./config/PhysiCell_settings.xml",
            "cell_definitions/cell_definition/phenotype/secretion/substrate/secretion_rate",
            name_cell_def="epi_inter",
            substrate="div_inhib",
        )
    )

    """print(modify_cell_data(
    file_path='./addons/start_and_stop/start_and_stop_scripts/cell_data.txt', 
    condition_key='type', 
    condition_val=1, 
    target_path='Secretion/Secretion 1/Secretion_Rate', 
    new_value=9.99
))"""
