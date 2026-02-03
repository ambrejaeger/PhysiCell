import re
import os

def modify_cell_data(file_path, condition_key, condition_val, target_path, new_value):
    """
    Modifies cell data in a custom text format based on a condition.
    
    :param file_path: Path to the .txt file
    :param condition_key: The key to check
    :param condition_val: The value the key must match
    :param target_path: List or string path
    :param new_value: The value to write
    """
    if not os.path.exists(file_path):
        return False

    with open(file_path, 'r') as f:
        content = f.read()

    # Split the file into individual Cell blocks
    # We use a lookahead to split by 'Cell' but keep the word 'Cell'
    blocks = re.split(r'(^Cell\n)', content, flags=re.MULTILINE)
    
    # Reassemble blocks (split creates empty entries and separates 'Cell' from its body)
    refined_blocks = []
    for i in range(1, len(blocks), 2):
        refined_blocks.append(blocks[i] + blocks[i+1])

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
        # Target path format: "Section/SubSection/Key"
        path_parts = target_path.split('/')
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
            if i > search_range_start and ":" in lines[i] and not any(s in lines[i] for s in [target_key] + sections):
                # This is a simple heuristic to ensure we don't jump into the next block section
                pass 
            
            if f"{target_key}:" in lines[i]:
                lines[i] = f"{target_key}: {new_value}"
                modifications_made = True
                break
        
        refined_blocks[b_idx] = "\n".join(lines) + "\n"

    # Save the file if changes were made
    if modifications_made:
        with open(file_path, 'w') as f:
            f.write("".join(refined_blocks))
        return True
    
    return False   

if __name__ == "__main__":
    print(modify_cell_data(
    file_path='./addons/start_and_stop/start_and_stop_scripts/cell_data.txt', 
    condition_key='type', 
    condition_val=1, 
    target_path='Secretion/Secretion 1/Secretion_Rate', 
    new_value=9.99
))