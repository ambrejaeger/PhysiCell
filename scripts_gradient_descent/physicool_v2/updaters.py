"""A module to create model updater functions for the PhysiCOOL black-box."""
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, Union, Callable, Optional
from dataclasses import dataclass, field

import physicool_v2.datatypes as dt
from physicool_v2.config import ConfigFileParser, CellRuleFileParser

def update_all(ruleset: dt.RuleSet, cell_data: dt.CellParameters, new_rules_values: list[tuple[int, str, float]], new_cell_values: Dict[str, float], death_model_name: str = "apoptosis"):
    update_ruleset(ruleset, new_rules_values)
    update_cell_values(cell_data, new_cell_values, death_model_name)

def update_ruleset(ruleset: dt.RuleSet, new_values: list[tuple[int, str, float]]):
    """Updates the numerical values for a RuleSet class."""
    for value in new_values:
        rule_index, parameter, new_value = value
        print(len(ruleset.rules))
        if rule_index < len(ruleset.rules):
            rule = ruleset.rules[rule_index]
            if hasattr(rule, parameter):
                setattr(rule, parameter, new_value)
            else:
                raise ValueError(f"Parameter {parameter} not found in the Rule class.")
        else:
            raise ValueError(f"Rule index {rule_index} is out of range for the ruleset.")
        
    return

CellUpdaterFunction = Callable[[dt.CellParameters, Dict[str, float]], None]

def update_cell_values(cell_data: dt.CellParameters, new_values: Dict[str, float], death_model_name: str = "apoptosis"):
    """
    Updates the numerical values for parameters of all type.

    Parameters
    ----------
    cell_data:
        The cell data structure to be modified.
    new_values:
        The new values to be written to the CellParameters class. 
    """
    update_cycle_values(cell_data, new_values)
    update_volume_values(cell_data, new_values)
    update_motility_values(cell_data, new_values)
    update_mechanics_values(cell_data, new_values)
    update_death_values(cell_data, new_values, death_model_name)


def update_cycle_values(cell_data: dt.CellParameters, new_values: Dict[str, float]):
    """
    Updates the numerical values for the Cycle class.

    Parameters
    ----------
    cell_data:
        The cell data structure to be modified.
    new_values:
        The new values to be written to the Cycle section
        of the CellParameters class. Keys should follow the pattern "phase_{i}"
        where i represents the id of the phase to be updated. All transition
        rates/durations must be specified.
    """
    if cell_data.cycle.phase_durations:
        if len(cell_data.cycle.phase_durations) != len(new_values):
            raise ValueError(
                "The passed values do not match the number of rates/durations."
            )

        cell_data.cycle.phase_durations = [
            new_values[f"phase_{i}"] for i, _ in enumerate(new_values)
        ]

    if cell_data.cycle.phase_transition_rates:
        if len(cell_data.cycle.phase_transition_rates) != len(new_values):
            raise ValueError(
                "The passed values do not match the number of rates/durations."
            )

        cell_data.cycle.phase_transition_rates = [
            new_values[f"phase_{i}"] for i, _ in enumerate(new_values)
        ]


def update_volume_values(cell_data: dt.CellParameters, new_values: Dict[str, float]):
    """
    Updates the numerical values for the Volume class.

    Parameters
    ----------
    cell_data:
        The cell data structure to be modified.
    new_values:
        The new values to be written to the volume section
        of the CellParameters class. Keys should be the same as those in the XML file,
        but it is not required to include all the keys.
    """
    if "total" in new_values.keys():
        cell_data.volume.total = new_values["total"]
    if "fluid_fraction" in new_values.keys():
        cell_data.volume.fluid_fraction = new_values["fluid_fraction"]
    if "nuclear" in new_values.keys():
        cell_data.volume.nuclear = new_values["nuclear"]
    if "fluid_change_rate" in new_values.keys():
        cell_data.volume.nuclear = new_values["fluid_change_rate"]
    if "cytoplasmic_biomass_change_rate" in new_values.keys():
        cell_data.volume.cytoplasmic_biomass_change_rate = new_values[
            "cytoplasmic_biomass_change_rate"
        ]
    if "nuclear_biomass_change_rate" in new_values.keys():
        cell_data.volume.nuclear_biomass_change_rate = new_values[
            "nuclear_biomass_change_rate"
        ]
    if "calcified_fraction" in new_values.keys():
        cell_data.volume.calcified_fraction = new_values["calcified_fraction"]
    if "calcification_rate" in new_values.keys():
        cell_data.volume.calcification_rate = new_values["calcification_rate"]
    if "relative_rupture_volume" in new_values.keys():
        cell_data.volume.relative_rupture_volume = new_values["relative_rupture_volume"]


def update_motility_values(cell_data: dt.CellParameters, new_values: Dict[str, float]):
    """
    Updates the numerical values for the Motility class.

    Parameters
    ----------
    cell_data:
        The cell data structure to be modified.
    new_values:
        The new values to be written to the motility section
        of the CellParameters class. Keys should be the same as those in the XML file,
        but it is not required to include all the keys.
    """
    if "speed" in new_values.keys():
        cell_data.motility.speed = new_values["speed"]
    if "persistence_time" in new_values.keys():
        cell_data.motility.persistence_time = new_values["persistence_time"]
    if "migration_bias" in new_values.keys():
        cell_data.motility.migration_bias = new_values["migration_bias"]


def update_mechanics_values(cell_data: dt.CellParameters, new_values: Dict[str, float]):
    """
    Updates the numerical values for the Mechanics class.

    Parameters
    ----------
    cell_data:
        The cell data structure to be modified.
    new_values:
        The new values to be written to the mechanics section
        of the CellParameters class. Keys should be the same as those in the XML file,
        but it is not required to include all the keys.
    """
    if "cell_cell_adhesion_strength" in new_values.keys():
        cell_data.mechanics.cell_cell_adhesion_strength = new_values[
            "cell_cell_adhesion_strength"
        ]
    if "cell_cell_repulsion_strength" in new_values.keys():
        cell_data.mechanics.cell_cell_repulsion_strength = new_values[
            "cell_cell_repulsion_strength"
        ]
    if "relative_maximum_adhesion_distance" in new_values.keys():
        cell_data.mechanics.relative_maximum_adhesion_distance = new_values[
            "relative_maximum_adhesion_distance"
        ]

def update_death_values(cell_data: dt.CellParameters, new_values: Dict[str, float], death_model_name: str = "apoptosis"):
    """
    Updates the numerical values for the death section of the CellParameters class.

    Parameters
    ----------
    cell_data:
        The cell data structure to be modified.
    
    new_values:
        The new values to be written to the death section
        of the CellParameters class. Keys should be the same as those in the XML file,
        but it is not required to include all the keys.
    """
    model_index = -1
    for i in range(len(cell_data.death)):
        if cell_data.death[i].name.lower() == death_model_name.lower():
            model_index = i
    if model_index == -1:
        raise ValueError(f"Death model {death_model_name} not found in the cell data.")
    
    if "death_rate" in new_values.keys():
        cell_data.death[model_index].death_rate = new_values["death_rate"]
    if "unlysed_fluid_change_rate" in new_values.keys():
        cell_data.death[model_index].unlysed_fluid_change_rate = new_values[
            "unlysed_fluid_change_rate"
        ]
    if "lysed_fluid_change_rate" in new_values.keys():
        cell_data.death[model_index].lysed_fluid_change_rate = new_values[
            "lysed_fluid_change_rate"
        ]
    if "cytoplasmic_biomass_change_rate" in new_values.keys():
        cell_data.death[model_index].cytoplasmic_biomass_change_rate = new_values[
            "cytoplasmic_biomass_change_rate"
        ]
    if "nuclear_biomass_change_rate" in new_values.keys():
        cell_data.death[model_index].nuclear_biomass_change_rate = new_values[
            "nuclear_biomass_change_rate"
        ]
    if "calcification_rate" in new_values.keys():
        cell_data.death[model_index].calcification_rate = new_values["calcification_rate"]
    if "relative_rupture_volume" in new_values.keys():
        cell_data.death[model_index].relative_rupture_volume = new_values["relative_rupture_volume"]

    for key in new_values.keys():
        if key.startswith("phase_durations"):
            index = key.split('[')[1].split(']')[0]
            if cell_data.death[model_index].phase_durations:
                cell_data.death[model_index].phase_durations[int(index)] = new_values[key]
            else:
                raise ValueError("Phase durations not defined for this death model in the cell data.")
        if key.startswith("phase_transition_rates"):
            index = key.split('[')[1].split(']')[0]
            if cell_data.death[model_index].phase_transition_rates:
                cell_data.death[model_index].phase_transition_rates[int(index)] = new_values[key]
            else:
                raise ValueError("Phase transition rates not defined for this death model in the cell data.")
    return
    
@dataclass
class ParamsUpdater(ABC):
    config_path: Union[str, Path]
    parser: ConfigFileParser = field(init=False)

    def __post_init__(self):
        """Creates the ConfigFileParser instance to be accessed by the class."""
        self.parser = ConfigFileParser(path=self.config_path)

    @abstractmethod
    def update(self, new_values: Dict[str, float]) -> None:
        """Updates the XML file with the values passed as input."""
        pass


@dataclass
class CellUpdater(ParamsUpdater):
    updater_function: CellUpdaterFunction
    cell_definition_name: str = "default"

    def update(self, new_values: Dict[str, float]) -> None:
        """Updates the XML file with the values passed as input."""
        cell_data = self.parser.read_cell_data(name=self.cell_definition_name)
        self.updater_function(cell_data, new_values)
        self.parser.write_cell_params(cell_data=cell_data)


def update_substance_values(substance: dt.Substance, new_values: Dict[str, float]):
    """
    Updates the numerical values for a Substance class (microenvironment).

    Parameters
    ----------
    substance
        The substance data to be updated.
    new_values
        The new values to be written to the substance class. Keys should be the same
        as those in the XML file, but it is not required to include all the keys.
    """
    if "diffusion_coefficient" in new_values.keys():
        substance.diffusion_coefficient = new_values["diffusion_coefficient"]
    if "decay_rate" in new_values.keys():
        substance.decay_rate = new_values["decay_rate"]
    if "initial_condition" in new_values.keys():
        substance.initial_condition = new_values["initial_condition"]
    if "dirichlet_boundary_condition" in new_values.keys():
        substance.dirichlet_boundary_condition = new_values[
            "dirichlet_boundary_condition"
        ]


@dataclass
class MicroenvironmentUpdater(ParamsUpdater):
    substance_name: str

    def update(self, new_values: Dict[str, float]) -> None:
        """Updates the XML file with the values passed as input."""
        substances = self.parser.read_me_params()
        substance = [
            substance
            for substance in substances
            if substance.name == self.substance_name
        ][0]
        update_substance_values(substance=substance, new_values=new_values)
        self.parser.write_substance_params(substance)

@dataclass
class RulesetUpdater:
    ruleset_parser: CellRuleFileParser = field(init=False)
    cfp: ConfigFileParser

    def __post_init__(self):
        """Creates the RuleFileParser instance to be accessed by the class."""
        ruleset_filepath = Path.cwd() / Path(self.cfp.read_ruleset_params().folder) / self.cfp.read_ruleset_params().filename       
        self.ruleset_parser = CellRuleFileParser(path=ruleset_filepath)

    def update(self, new_values: list[tuple[int, str, float]]) -> None:
        """Updates the ruleset file with the values passed as input."""
        ruleset = self.ruleset_parser.read_ruleset()
        update_ruleset(ruleset, new_values)
        # write back each rule to the cell rules file
        for i, rule in enumerate(ruleset.rules):
            self.ruleset_parser.write_cell_rules(i, rule)