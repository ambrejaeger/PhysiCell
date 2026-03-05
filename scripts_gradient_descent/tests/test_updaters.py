import unittest

import physicool_v2.datatypes as dt
from physicool_v2 import updaters
from physicool_v2 import config
from tests.configdata import *
import shutil

CELL_DATA = {
    "name": "default",
    "cycle": {
        "code": 6.0,
        "phase_durations": [300.0, 480.0, 240.0, 60.0],
        "phase_transition_rates": None,
    },
    "death": [
        {
            "code": 100.0,
            "name": "apoptosis",
            "death_rate": 5.31667e-05,
            "phase_durations": [516.0],
            "phase_transition_rates": None,
            "unlysed_fluid_change_rate": 0.05,
            "lysed_fluid_change_rate": 0.0,
            "cytoplasmic_biomass_change_rate": 0.0166667,
            "nuclear_biomass_change_rate": 0.00583333,
            "calcification_rate": 0.0,
            "relative_rupture_volume": 2.0,
        },
        {
            "code": 101.0,
            "name": "necrosis",
            "death_rate": 0.0,
            "phase_durations": [0.0, 86400.0],
            "phase_transition_rates": None,
            "unlysed_fluid_change_rate": 0.05,
            "lysed_fluid_change_rate": 0.0,
            "cytoplasmic_biomass_change_rate": 0.0166667,
            "nuclear_biomass_change_rate": 0.00583333,
            "calcification_rate": 0.0,
            "relative_rupture_volume": 2.0,
        },
    ],
    "volume": {
        "total": 2494.0,
        "fluid_fraction": 0.75,
        "nuclear": 540.0,
        "fluid_change_rate": 0.05,
        "cytoplasmic_biomass_change_rate": 0.0045,
        "nuclear_biomass_change_rate": 0.0055,
        "calcified_fraction": 0.0,
        "calcification_rate": 0.0,
        "relative_rupture_volume": 2.0,
    },
    "mechanics": {
        "cell_cell_adhesion_strength": 0.4,
        "cell_cell_repulsion_strength": 10.0,
        "relative_maximum_adhesion_distance": 1.25,
        "set_relative_equilibrium_distance": 1.8,
        "set_absolute_equilibrium_distance": 15.12,
    },
    "motility": {
        "speed": 1.0,
        "persistence_time": 1.0,
        "migration_bias": 0.5,
        "motility_enabled": False,
        "use_2d": True,
        "chemotaxis_enabled": False,
        "chemotaxis_substrate": "substrate",
        "chemotaxis_direction": 1.0,
    },
    "secretion": [
        {
            "name": "substrate",
            "secretion_rate": 0.0,
            "secretion_target": 1.0,
            "uptake_rate": 0.0,
            "net_export_rate": 0.0,
        },
        {
            "name": "oxygen",
            "secretion_rate": 0.0,
            "secretion_target": 1.0,
            "uptake_rate": 0.0,
            "net_export_rate": 0.0,
        },
    ],
    "custom": [{"name": "sample", "value": 1.0}],
}


EXPECTED_DEATH = [dt.Death(
            code=100.0,
            name="apoptosis",
            death_rate=0.1,
            phase_durations=[500.0],
            phase_transition_rates= None,
            unlysed_fluid_change_rate=0.1,
            lysed_fluid_change_rate=0.1,
            cytoplasmic_biomass_change_rate=0.1,
            nuclear_biomass_change_rate=0.1,
            calcification_rate=0.1,
            relative_rupture_volume=0.1,
        ),
        dt.Death(
            code=101.0,
            name="necrosis",
            death_rate=0.0,
            phase_durations=[0.0, 86400.0],
            phase_transition_rates=None,
            unlysed_fluid_change_rate=0.05,
            lysed_fluid_change_rate=0.0,
            cytoplasmic_biomass_change_rate=0.0166667,
            nuclear_biomass_change_rate=0.00583333,
            calcification_rate=0.0,
            relative_rupture_volume=2.0,
        )]

EXPECTED_CYCLE = dt.Cycle(
    code=6.0, phase_durations=[20.0, 180.0, 240.0, 60.0], phase_transition_rates=None
)

EXPECTED_VOLUME = dt.Volume(
    total=2494.0,
    fluid_fraction=0.75,
    nuclear=540.0,
    fluid_change_rate=0.05,
    cytoplasmic_biomass_change_rate=0.0045,
    nuclear_biomass_change_rate=0.0055,
    calcified_fraction=0.0,
    calcification_rate=0.0,
    relative_rupture_volume=2.0,
)

EXPECTED_MOTILITY = dt.Motility(
    speed=5.0,
    persistence_time=10.0,
    migration_bias=1.0,
    motility_enabled=False,
    use_2d=True,
    chemotaxis_enabled=False,
    chemotaxis_substrate="substrate",
    chemotaxis_direction=1.0,
)

EXPECTED_MOTILITY_2 = EXPECTED_MOTILITY.copy(deep=True)
EXPECTED_MOTILITY_2.migration_bias = 0.5

EXPECTED_DATA = dt.CellParameters(
    name="default",
    cycle=EXPECTED_CYCLE,
    death=EXPECTED_DEATH,
    volume=EXPECTED_VOLUME,
    mechanics=dt.Mechanics(**EXPECTED_MECHANICS_READ),
    motility=EXPECTED_MOTILITY,
    secretion=[dt.Secretion(**EXPECTED_SECRETION_READ_SUBSTRATE)],
    custom=[dt.CustomData(**EXPECTED_CUSTOM_READ[0])]
)

EXPECTED_RULESET = dt.RuleSet(rules=[
    dt.Rule(
        type='epi_basal',
        signal='contact with membrane',
        monotony='decreases',
        behavior='transform to epi_inter',
        base_value=float(0.0),
        half_max=float(0.1),
        hill_power=float(4),
        apply_to_dead=bool(0)
    ),
    dt.Rule(
        type='epi_basal',
        signal='div_inhib',
        monotony='decreases',
        behavior='exit from cycle phase 0',
        base_value=float(0.0),
        half_max=float(1.0),
        hill_power=float(4),
        apply_to_dead=bool(0)
    ),
    dt.Rule(
        type='epi_inter',
        signal='contact with epi_inter',
        monotony='decreases',
        behavior='apoptosis',
        base_value=float(1.0),
        half_max=float(0.001),
        hill_power=float(4),
        apply_to_dead=bool(1)
    )
])

class UpdaterFunctionsTest(unittest.TestCase):
    def test_update_all(self):
        """Asserts that all the parameters are correctly updated with update all."""
        data = dt.CellParameters(**CELL_DATA)
        new_values = {
            "phase_0": 20.0,
            "phase_1": 180.0,
            "phase_2": 240.0,
            "phase_3": 60.0,
            "phase_durations[0]": 500.0,
            "death_rate": 0.1,
            "unlysed_fluid_change_rate": 0.1,
            "lysed_fluid_change_rate": 0.1,
            "cytoplasmic_biomass_change_rate": 0.1,
            "nuclear_biomass_change_rate": 0.1,
            "calcification_rate": 0.1,
            "relative_rupture_volume": 0.1,
            "speed": 5.0,
            "persistence_time": 10.0,
            "migration_bias": 1.0,
        }
        updaters.update_death_values(cell_data=data, new_values=new_values)
        self.assertEqual(EXPECTED_DATA, data)

    def test_cell_rule_updater_function(self):
        #write a EXPECTED_RULESET

        cfp = config.ConfigFileParser(CONFIG_PATH)
        rulefileupdater = updaters.RulesetUpdater(cfp)
        print('ruleset_path: ', rulefileupdater.ruleset_parser.cell_rule_file)
        ruleset_path = Path(rulefileupdater.ruleset_parser.cell_rule_file)

        # Make a copy of the file
        new_path = ruleset_path.parent / f"test_{ruleset_path.name}"
        shutil.copy2(ruleset_path, new_path)
        rulefileupdater.ruleset_parser.cell_rule_file = new_path

        #modify the copy of the cell rule file
        new_values = [(2, 'base_value', 1.0), (2, 'apply_to_dead', 1)]

        #write the new values to the new copy of the cell rule file
        rulefileupdater.update(new_values)

        #read the copy
        mod_ruleset = rulefileupdater.ruleset_parser.read_ruleset()
        print(mod_ruleset)
        self.assertEqual(EXPECTED_RULESET, mod_ruleset)

    def test_death_updater_function(self):
        """Asserts that the cycle parameters are correctly updated."""
        data = dt.CellParameters(**CELL_DATA)
        new_phase_values = {
            "phase_durations[0]": 500.0,
            "death_rate": 0.1,
            "unlysed_fluid_change_rate": 0.1,
            "lysed_fluid_change_rate": 0.1,
            "cytoplasmic_biomass_change_rate": 0.1,
            "nuclear_biomass_change_rate": 0.1,
            "calcification_rate": 0.1,
            "relative_rupture_volume": 0.1
        }
        updaters.update_death_values(cell_data=data, new_values=new_phase_values)
        self.assertEqual(EXPECTED_DEATH, data.death)

    def test_cycle_updater_function(self):
        """Asserts that the cycle parameters are correctly updated."""
        data = dt.CellParameters(**CELL_DATA)
        new_cycle_values = {
            "phase_0": 20.0,
            "phase_1": 180.0,
            "phase_2": 240.0,
            "phase_3": 60.0,
        }
        updaters.update_cycle_values(cell_data=data, new_values=new_cycle_values)
        self.assertEqual(EXPECTED_CYCLE, data.cycle)

    def test_cycle_updater_function_wrong_length(self):
        """Asserts that the cycle parameters are correctly updated."""
        data = dt.CellParameters(**CELL_DATA)
        new_cycle_values = {
            "phase_0": 20.0,
            "phase_1": 180.0,
            "phase_2": 240.0,
        }
        self.assertRaises(
            ValueError, updaters.update_cycle_values, data, new_values=new_cycle_values
        )

    def test_motility_updater_function(self):
        """Asserts that the motility parameters are correctly updated."""
        data = dt.CellParameters(**CELL_DATA)
        new_motility_values = {
            "speed": 5.0,
            "persistence_time": 10.0,
            "migration_bias": 1.0,
        }
        updaters.update_motility_values(cell_data=data, new_values=new_motility_values)
        self.assertEqual(EXPECTED_MOTILITY, data.motility)

    def test_motility_updater_function_incomplete(self):
        """Asserts that the motility parameters are correctly updated when not all parameters are defined."""
        data = dt.CellParameters(**CELL_DATA)
        new_motility_values = {"speed": 5.0, "persistence_time": 10.0}
        updaters.update_motility_values(cell_data=data, new_values=new_motility_values)
        self.assertEqual(EXPECTED_MOTILITY_2, data.motility)


if __name__ == "__main__":
    unittest.main()
