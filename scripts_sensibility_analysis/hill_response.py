#A short script to shake the aspect of the hill function in the rules file and see 
#how it changes the response curve.

import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET


def down_hill_function(x,b0,bm,hfm,hp):
    """
    b0: value in the absence of signal, it is read from the PhysiCell settings file
    bm: minimum value when there is a signal, it is read from the rule file
    """
    x = np.asarray(x)
    y = b0 + (bm - b0) * (x/hfm)**hp / (1 + (x/hfm)**hp)
    return x, y

def up_hill_function(x,b0,bm,hfm,hp):
    """
    b0: value in the absence of signal, it is read from the PhysiCell settings file
    bm: maximum value when there is a signal, it is read from the rule file
    """
    x = np.asarray(x)
    y = b0 + (bm - b0) * (x/hfm)**hp / (1 + (x/hfm)**hp)
    return x, y

def rule_to_hill_function(rule, x, b0):
    if len(rule) != 7: 
        rule_type = rule[2]  # 3rd element
        bm = float(rule[4])  # 5th element
        halfmax = float(rule[5])  # 6th element
        hillpower = float(rule[6])  # 7th element
    
        # Select appropriate function
        if rule_type == 'decreases':
            if b0 < bm:
                raise ValueError(f"b0 ({b0}) should be greater than bm ({bm}) for a decreasing rule.")
            return down_hill_function(x,b0,bm,halfmax,hillpower)
        elif rule_type == 'increases':
            if b0 > bm:
                raise ValueError(f"b0 ({b0}) should be less than bm ({bm}) for an increasing rule.")
            return up_hill_function(x,b0,bm,halfmax, hillpower)
    else:
        raise ValueError("Rule is not of the correct size")
    
def read_rules_file(filepath):
    rules = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('//'):
                continue

            rule_parts = line.split(',')
            rules.append(rule_parts)
    
    return rules


def read_value_in_settings_file(setting_file,path):
    """setting_file: path to the setting file
       path: path to the variable of interest in the xml setting file
       it should be in the form path/to/the/value and attribute are accessed 
       with [@attribute='value'], for instance:
       path = "cell_definitions/cell_definition[@name='epi_basal']/phenotype/death/model[@name='apoptosis']/death_rate"
    """
    tree = ET.parse(setting_file)
    root = tree.getroot()
    value = root.find(path).text
    return value



if __name__ == "__main__":
    path = "/home/ajaeger/Documents/PhysiCell/config/cell_rules.csv"
    setting_file = "/home/ajaeger/Documents/PhysiCell/config/PhysiCell_settings.xml"
    rules = read_rules_file(path)

    print(rules)

    plt.figure(figsize=(8, 5))
    x = np.linspace(0, 1, 100)
    #path_to_value = "cell_definitions/cell_definition[@name='epi_inter']/phenotype/death/model[@name='apoptosis']/death_rate"
    path_to_value = "cell_definitions/cell_definition[@name='epi_basal']/phenotype/cycle/phase_transition_rates/rate"
    
    b0 = read_value_in_settings_file(setting_file,path_to_value)
    print("bO is equal to: ", b0)
    x, y = rule_to_hill_function(rules[1],x,float(b0))

    plt.plot(x, y, label=f'b={rules[1][4]}, hfm={rules[1][5]}, hp={rules[1][6]}', linewidth=2)
    rules[1][5] = 3
    rules[1][6] = 0.5
    x, y = rule_to_hill_function(rules[1],x,float(b0))
    plt.plot(x, y, label=f'b={1}, hfm={rules[1][5]}, hp={rules[1][6]}', linewidth=2)
    plt.show()