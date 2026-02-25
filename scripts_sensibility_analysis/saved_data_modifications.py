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
        print(f"Error: File '{file_path}' does not exist.")
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
                    break

        if not cell_matches:
            continue

        # 2. If matched, navigate the target path
        # Target path format: "Section/SubSection/Key" optionally with an
        # index for list entries, e.g. "Section/Key[2]" to modify the 3rd
        # whitespace-separated token on that line.
        path_parts = target_path.split("/")
        raw_target_key = path_parts[-1]
        sections = path_parts[:-1]

        # Detect optional index in the target key (e.g. Key[2])
        list_index = None
        m = re.match(r"^(?P<key>.+)\[(?P<idx>\d+)\]$", raw_target_key)
        if m:
            target_key = m.group("key")
            list_index = int(m.group("idx"))
        else:
            target_key = raw_target_key

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
                # Extract existing value part after the colon
                parts = lines[i].split(":", 1)
                if len(parts) == 2:
                    prefix = parts[0].strip()
                    value_part = parts[1].strip()
                else:
                    prefix = parts[0].strip()
                    value_part = ""

                # If an index was requested and the value is a whitespace-separated list,
                # replace only that token and keep others unchanged.
                if list_index is not None and value_part != "":
                    tokens = value_part.split()
                    if 0 <= list_index < len(tokens):
                        tokens[list_index] = str(new_value)
                        new_value_str = " ".join(tokens)
                    else:
                        # If index out of range, pad with zeros up to that index
                        # then set the requested index.
                        while len(tokens) <= list_index:
                            tokens.append("0")
                        tokens[list_index] = str(new_value)
                        new_value_str = " ".join(tokens)
                else:
                    # Replace the whole value as before
                    new_value_str = str(new_value)

                lines[i] = f"{prefix}: {new_value_str}"
                modifications_made = True
                break

        refined_blocks[b_idx] = "\n".join(lines) + "\n"

    # Save the file if changes were made
    if modifications_made:
        with open(file_path, "w") as f:
            f.write("".join(refined_blocks))
        return True

    return False


if __name__ == "__main__":
 
    print(modify_cell_data(
    file_path='./addons/start_and_stop/start_and_stop_scripts/cell_data.txt', 
    condition_key='type', 
    condition_val="1", 
    target_path='Mechanics/cell_adhesion_affinities[1]', 
    new_value="9.99"
))
