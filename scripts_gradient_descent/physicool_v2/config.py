"""
This module enables users to modify their PhysiCell XML config file
and to modify the Cell rule CSV file associated to it
"""
from pathlib import Path
from xml.etree import ElementTree
from typing import List, Union

import physicool_v2.datatypes as dt
from physicool_v2 import pcxml
import csv

class ConfigFileParser:
    """
    A class that acts as an interface between the user and the XML config file.

    Parameters
    ----------
    path
        The path to the configuration file to be read/written by the parser.
    """

    def __init__(
        self, path: Union[str, Path] = Path("config/PhysiCell_settings.xml")
    ) -> None:
        if isinstance(path, str):
            path = Path(path)

        self.config_file = path
        self.tree = ElementTree.parse(path)

    def __repr__(self):
        return f"ConfigFileParser(config_file={self.config_file})"

    @property
    def cell_definitions_list(self) -> List[str]:
        """Returns a list with the names of the cell definitions in the XML file."""
        root = self.tree.getroot()
        cell_definitions = root.find("cell_definitions").findall("cell_definition")

        return [definition.attrib["name"] for definition in cell_definitions]

    @property
    def me_substance_list(self) -> List[str]:
        """Returns a list with the names of the microenvironment substances in the XML file."""
        root = self.tree.getroot()
        substances = root.find("microenvironment/domain/variables").findall("variable")

        return [substance.attrib["name"] for substance in substances]
    
    def read_ruleset_params(self) -> dt.RuleSetParams:
        """Returns the <cell_rules> data from the XML file."""
        return dt.RuleSetParams(**pcxml.parse_ruleset_params(tree=self.tree, path="cell_rules/rulesets/ruleset"))
    
    def read_domain_params(self) -> dt.Domain:
        """Returns the <domain> data from the XML file."""
        return dt.Domain(**pcxml.parse_domain(tree=self.tree, path="domain"))

    def read_overall_params(self) -> dt.Overall:
        """Returns the <overall> data from the XML file."""
        return dt.Overall(**pcxml.parse_overall(tree=self.tree, path="overall"))

    def read_me_params(self) -> List[dt.Substance]:
        """Returns the <microenvironment_setup> data form the XML file."""
        return [
            dt.Substance(**substance)
            for substance in pcxml.parse_microenvironment(
                tree=self.tree, path="microenvironment_setup"
            )
        ]

    def read_cycle_params(self, name: str) -> dt.Cycle:
        """
        Returns the <cycle> data for a given cell definition from the XML file.

        Parameters
        ----------
        name
            A string with the name of the cell definition to be read
        """
        stem = f"cell_definitions/cell_definition[@name='{name}']/phenotype/cycle"
        return dt.Cycle(**pcxml.parse_cycle(self.tree, path=stem))

    def read_death_params(self, name: str) -> List[dt.Death]:
        """
        Returns the <death> data for a given cell definition from the XML file.

        Parameters
        ----------
        name
            A string with the name of the cell definition to be read
        """
        stem = f"cell_definitions/cell_definition[@name='{name}']/phenotype/death"
        return [dt.Death(**model) for model in pcxml.parse_death(self.tree, path=stem)]

    def read_volume_params(self, name: str) -> dt.Volume:
        """
        Returns the <volume> data for a given cell definition from the XML file.

        Parameters
        ----------
        name
            A string with the name of the cell definition to be read
        """
        stem = f"cell_definitions/cell_definition[@name='{name}']/phenotype/volume"
        return dt.Volume(**pcxml.parse_volume(tree=self.tree, path=stem))

    def read_mechanics_params(self, name: str) -> dt.Mechanics:
        """
        Returns the <mechanics> data for a given cell definition from the XML file.

        Parameters
        ----------
        name
            A string with the name of the cell definition to be read
        """
        stem = f"cell_definitions/cell_definition[@name='{name}']/phenotype/mechanics"
        return dt.Mechanics(**pcxml.parse_mechanics(tree=self.tree, path=stem))

    def read_motility_params(self, name: str) -> dt.Motility:
        """
        Returns the <motility> data for a given cell definition from the XML file.

        Parameters
        ----------
        name
            A string with the name of the cell definition to be read
        """
        stem = f"cell_definitions/cell_definition[@name='{name}']/phenotype/motility"
        return dt.Motility(**pcxml.parse_motility(tree=self.tree, path=stem))

    def read_secretion_params(self, name: str) -> List[dt.Secretion]:
        """
        Returns the <secretion> data for a given cell definition from the XML file.

        Parameters
        ----------
        name
            A string with the name of the cell definition to be read
        """
        stem = f"cell_definitions/cell_definition[@name='{name}']/phenotype/secretion"
        return [
            dt.Secretion(**substance)
            for substance in pcxml.parse_secretion(tree=self.tree, path=stem)
        ]

    def read_custom_data(self, name: str) -> List[dt.CustomData]:
        """
        Returns the <custom_data> data for a given cell definition from the XML file.

        Parameters
        ----------
        name
            A string with the name of the cell definition to be read
        """
        stem = f"cell_definitions/cell_definition[@name='{name}']/custom_data"
        return [
            dt.CustomData(**custom) for custom in pcxml.parse_custom(self.tree, stem)
        ]

    def read_cell_data(self, name: str = "default") -> dt.CellParameters:
        """
        Reads all the fields for a given cell definition into a custom cell data type.

        Parameters
        ----------
        name
            The name of the cell definition to be read.

        Returns
        -------
        dt.CellParameters
            A custom cell data type that contains all the parameters of a cell definition.

        Raises
        ------
        ValueError
            If the passed cell definition is not defined in the config file.
        """
        if name not in self.cell_definitions_list:
            raise ValueError("Invalid cell definition")

        # Read and save the cell data
        cycle = self.read_cycle_params(name)
        death = self.read_death_params(name)
        volume = self.read_volume_params(name)
        mechanics = self.read_mechanics_params(name)
        motility = self.read_motility_params(name)
        secretion = self.read_secretion_params(name)
        custom = self.read_custom_data(name)

        return dt.CellParameters(
            name=name,
            cycle=cycle,
            death=death,
            volume=volume,
            mechanics=mechanics,
            motility=motility,
            secretion=secretion,
            custom=custom,
        )

    def read_user_params(self):
        """Returns the <user_parameters> data  from the XML file."""
        return [
            dt.CustomData(**custom)
            for custom in pcxml.parse_custom(self.tree, "user_parameters")
        ]

    def write_domain_params(self, domain: dt.Domain, update_file: bool = True) -> None:
        """
        Writes the passed <domain> data to the XML tree and file.

        Parameters
        ----------
        domain
            The new data values to be written.
        update_file
            If the values should be written to the file. If False, the values
            will only be changed in the XML tree.
        """
        pcxml.write_domain(new_values=domain.dict(), tree=self.tree, path="domain")
        if update_file:
            self.tree.write(self.config_file)

    def write_overall_params(
        self, overall: dt.Overall, update_file: bool = True
    ) -> None:
        """
        Writes the passed <overall> data to the XML tree and file.

        Parameters
        ----------
        overall
            The new data values to be written.
        update_file
            If the values should be written to the file. If False, the values
            will only be changed in the XML tree.
        """
        pcxml.write_overall(new_values=overall.dict(), tree=self.tree, path="overall")
        if update_file:
            self.tree.write(self.config_file)

    def write_substance_params(
        self, substance: dt.Substance, update_file: bool = True
    ) -> None:
        """
        Writes the passed <microenvironment_setup> data for a given
        substance to the XML tree and file.

        Parameters
        ----------
        substance
            The new data values to be written.
        update_file
            If the values should be written to the file. If False, the values
            will only be changed in the XML tree.
        """
        pcxml.write_substance(
            new_values=substance.dict(),
            tree=self.tree,
            path="microenvironment_setup",
            name=substance.name,
        )
        if update_file:
            self.tree.write(self.config_file)

    def write_cycle_params(
        self, name: str, cycle: dt.Cycle, update_file: bool = True
    ) -> None:
        """
        Writes the passed <cycle> data to the XML tree and file.

        Parameters
        ----------
        name
            The name of the cell definition to be updated.
        cycle
            The new data values to be written.
        update_file
            If the values should be written to the file. If False, the values
            will only be changed in the XML tree.

        Raises
        -----
        ValueError
            When the passed cell definition name does not match any of the
            cell definitions in the file.
        """
        if name not in self.cell_definitions_list:
            raise ValueError("The passed cell definition is not in the XML file.")

        stem = f"cell_definitions/cell_definition[@name='{name}']/phenotype/cycle"
        pcxml.write_cycle(new_values=cycle.dict(), tree=self.tree, path=stem)
        if update_file:
            self.tree.write(self.config_file)

    def write_death_model_params(
        self, name: str, death: dt.Death, update_file: bool = True
    ) -> None:
        """
        Writes the passed <death> data to the XML tree and file.

        Parameters
        ----------
        name
            The name of the cell definition to be updated.
        death
            The new data values to be written.
        update_file
            If the values should be written to the file. If False, the values
            will only be changed in the XML tree.

        Raises
        -----
        ValueError
            When the passed cell definition name does not match any of the
            cell definitions in the file.
        """
        if name not in self.cell_definitions_list:
            raise ValueError("The passed cell definition is not in the XML file.")

        stem = f"cell_definitions/cell_definition[@name='{name}']/phenotype/death"
        pcxml.write_death_model(new_values=death.dict(), tree=self.tree, path=stem)
        if update_file:
            self.tree.write(self.config_file)

    def write_death_params(
        self, name: str, death: List[dt.Death], update_file: bool = True
    ) -> None:
        """
        Writes the passed <volume data to the XML tree and file.

        Parameters
        ----------
        name
            The name of the cell definition to be updated.
        death
            The new data values to be written.
        update_file
            If the values should be written to the file. If False, the values
            will only be changed in the XML tree.

        Raises
        -----
        ValueError
            When the passed cell definition name does not match any of the
            cell definitions in the file.
        """
        if name not in self.cell_definitions_list:
            raise ValueError("The passed cell definition is not in the XML file.")

        stem = f"cell_definitions/cell_definition[@name='{name}']/phenotype/death"
        for model in death:
            pcxml.write_death_model(new_values=model.dict(), tree=self.tree, path=stem)

        if update_file:
            self.tree.write(self.config_file)

    def write_volume_params(
        self, name: str, volume: dt.Volume, update_file: bool = True
    ) -> None:
        """
        Writes the passed <volume data to the XML tree and file.

        Parameters
        ----------
        name
            The name of the cell definition to be updated.
        volume
            The new data values to be written.
        update_file
            If the values should be written to the file. If False, the values
            will only be changed in the XML tree.

        Raises
        -----
        ValueError
            When the passed cell definition name does not match any of the
            cell definitions in the file.
        """
        if name not in self.cell_definitions_list:
            raise ValueError("The passed cell definition is not in the XML file.")

        stem = f"cell_definitions/cell_definition[@name='{name}']/phenotype/volume"
        pcxml.write_volume(new_values=volume.dict(), tree=self.tree, path=stem)
        if update_file:
            self.tree.write(self.config_file)

    def write_mechanics_params(
        self, name: str, mechanics: dt.Mechanics, update_file: bool = True
    ) -> None:
        """
        Writes the passed <mechanics> data to the XML tree and file.

        Parameters
        ----------
        name
            The name of the cell definition to be updated.
        mechanics
            The new data values to be written.
        update_file
            If the values should be written to the file. If False, the values
            will only be changed in the XML tree.

        Raises
        -----
        ValueError
            When the passed cell definition name does not match any of the
            cell definitions in the file.
        """
        if name not in self.cell_definitions_list:
            raise ValueError("The passed cell definition is not in the XML file.")

        stem = f"cell_definitions/cell_definition[@name='{name}']/phenotype/mechanics"
        pcxml.write_mechanics(new_values=mechanics.dict(), tree=self.tree, path=stem)
        if update_file:
            self.tree.write(self.config_file)

    def write_motility_params(
        self, name: str, motility: dt.Motility, update_file: bool = True
    ) -> None:
        """
        Writes the passed <motility> data to the XML tree and file.

        Parameters
        ----------
        name
            The name of the cell definition to be updated.
        motility
            The new data values to be written.
        update_file
            If the values should be written to the file. If False, the values
            will only be changed in the XML tree.

        Raises
        -----
        ValueError
            When the passed cell definition name does not match any of the
            cell definitions in the file.
        """
        if name not in self.cell_definitions_list:
            raise ValueError("The passed cell definition is not in the XML file.")

        stem = f"cell_definitions/cell_definition[@name='{name}']/phenotype/motility"
        pcxml.write_motility(new_values=motility.dict(), tree=self.tree, path=stem)
        if update_file:
            self.tree.write(self.config_file)

    def write_secretion_substance_params(
        self,
        name: str,
        secretion: dt.Secretion,
        substance: str,
        update_file: bool = True,
    ) -> None:
        """
        Writes the passed <secretion> data for a substance to the XML tree and file.

        Parameters
        ----------
        name
            The name of the cell definition to be updated.
        secretion
            The new data values to be written.
        substance
            The name of the secretion substance to be updated.
        update_file
            If the values should be written to the file. If False, the values
            will only be changed in the XML tree.

        Raises
        -----
        ValueError
            When the passed cell definition name does not match any of the
            cell definitions in the file.
        """
        if name not in self.cell_definitions_list:
            raise ValueError("The passed cell definition is not in the XML file.")

        stem = f"cell_definitions/cell_definition[@name='{name}']/phenotype/secretion"
        pcxml.write_secretion_substance(
            new_values=secretion.dict(), tree=self.tree, path=stem, name=substance
        )
        if update_file:
            self.tree.write(self.config_file)

    def write_secretion_params(
        self, name: str, secretion: List[dt.Secretion], update_file: bool = True
    ) -> None:
        """
        Writes the passed <secretion> data to the XML tree and file.

        Parameters
        ----------
        name
            The name of the cell definition to be updated.
        secretion
            The new data values to be written.
        update_file
            If the values should be written to the file. If False, the values
            will only be changed in the XML tree.

        Raises
        -----
        ValueError
            When the passed cell definition name does not match any of the
            cell definitions in the file.
        """
        if name not in self.cell_definitions_list:
            raise ValueError("The passed cell definition is not in the XML file.")

        stem = f"cell_definitions/cell_definition[@name='{name}']/phenotype/secretion"

        for substance in secretion:
            pcxml.write_secretion_substance(
                new_values=substance.model_dump(),
                tree=self.tree,
                path=stem,
                name=substance.name,
            )
        if update_file:
            self.tree.write(self.config_file)

    def write_custom_params(
        self, name: str, custom_data: List[dt.CustomData], update_file: bool = True
    ):
        """
        Writes the passed <custom_data> data to the XML tree and file.

        Parameters
        ----------
        name
            The name of the cell definition to be updated.
        custom_data
            The new data values to be written.
        update_file
            If the values should be written to the file. If False, the values
            will only be changed in the XML tree.

        Raises
        -----
        ValueError
            When the passed cell definition name does not match any of the
            cell definitions in the file.
        """
        if name not in self.cell_definitions_list:
            raise ValueError("The passed cell definition is not in the XML file.")

        stem = f"cell_definitions/cell_definition[@name='{name}']/custom_data"
        data = [variable.dict() for variable in custom_data]
        pcxml.write_custom_data(new_values=data, tree=self.tree, path=stem)
        if update_file:
            self.tree.write(self.config_file)

    def write_user_params(
        self, custom_data: List[dt.CustomData], update_file: bool = True
    ):
        """
        Writes the passed <user_parameters> data to the XML tree and file.

        Parameters
        ----------
        custom_data
            The new data values to be written.
        update_file
            If the values should be written to the file. If False, the values
            will only be changed in the XML tree.
        """
        data = [variable.dict() for variable in custom_data]
        pcxml.write_custom_data(new_values=data, tree=self.tree, path="user_parameters")
        if update_file:
            self.tree.write(self.config_file)

    def write_cell_params(self, cell_data: dt.CellParameters) -> None:
        """
        Writes the new parameters to the XML tree object and updates the XML file.

        Parameters
        ----------
        cell_data: dt.CellParameters
            The new cell parameters to be written to the XML file.
        """
        self.write_cycle_params(cell_data.name, cell_data.cycle, update_file=False)
        self.write_death_params(cell_data.name, cell_data.death, update_file=False)
        self.write_volume_params(cell_data.name, cell_data.volume, update_file=False)
        self.write_mechanics_params(
            cell_data.name, cell_data.mechanics, update_file=False
        )
        self.write_motility_params(
            cell_data.name, cell_data.motility, update_file=False
        )
        self.write_custom_params(cell_data.name, cell_data.custom, update_file=False)
        self.tree.write(self.config_file)

class CellRuleFileParser:
    """
    A class that acts as an interface between the user and the cell rules CSV file.

    Parameters
    ----------
    path
        The path to the cell rules file to be read/written by the parser.
    """

    def __init__(
        self, path: Union[str, Path] = Path("config/cell_rules.csv")
    ) -> None:
        if isinstance(path, str):
            path = Path(path)

        self.cell_rule_file = path
        self.ruleset = pcxml.parse_ruleset(path)


    def __repr__(self):
        return f"CellRuleFileParser(cell_rule_file={self.cell_rule_file})"
    
    def write(self, path: Union[str, Path]) -> None:
        """
        Write the ruleset to a CSV file, with one rule per line.
        
        Parameters
        ----------
        path
            The path where to write the file. If None, uses the original file path.
        """
        # Determine the output path
        output_path = Path(path) if path else self.cell_rule_file
        
        # Write each rule as a separate line
        with open(output_path, 'w', encoding='utf-8') as f:
            for rule in self.ruleset:
                if rule[7] == "True":
                    rule[7] = "1"
                elif rule[7] == "False":
                    rule[7] = "0"
                f.write(f"{','.join(rule)}\n")
    
    def read_ruleset(self) -> dt.RuleSet:
        """
        Parses a CSV file and returns a list of dictionaries with the specified schema.
        
        Schema order:
        type: str, signal: str, monotony: str, behavior: str, 
        base_value: float, half_max: float, hill_power: float, apply_to_dead: bool
        """
        rules = []
        
        with open(self.cell_rule_file, 'r') as file:
            reader = csv.reader(file)
            for i, row in enumerate(reader, 1):
                if not row or all(cell == '' for cell in row):  # Skip empty lines
                    continue
                    
                # Ensure we have all 8 columns
                if len(row) != 8:
                    print(f"Warning: Line {i} has {len(row)} columns, expected 8. Skipping.")
                    continue
                else:
                    rule = dt.Rule(**pcxml.parse_cell_rule(row))
                    print(rule)
                    rules.append(rule)
                    
        return dt.RuleSet(rules = rules)
    
    def write_cell_rules(self, rule_index: int, new_rule: dt.Rule, update_file: bool = True) -> None:
        print("The ruleset is: ", self.ruleset)
        if 0 <= rule_index < len(self.ruleset):
            print(new_rule)
            pcxml.write_cell_rule(rule_index, new_rule.model_dump(), self.ruleset)
        print("The ruleset is: ", self.ruleset)
        if update_file:
            self.write(self.cell_rule_file)

