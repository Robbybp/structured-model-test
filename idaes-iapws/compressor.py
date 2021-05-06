import pyomo.environ as pyo
from idaes.core import FlowsheetBlock

from idaes.generic_models.properties import iapws95
import idaes.generic_models.properties.helmholtz.helmholtz as hltz
from idaes.power_generation.unit_models.helm.compressor import (
        HelmIsentropicCompressor,
        )

from idaes.core.util.model_statistics import degrees_of_freedom


def make_model():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.prop_water = iapws95.Iapws95ParameterBlock(
            default={
                "phase_presentation": iapws95.PhaseType.MIX,
                "state_vars": iapws95.StateVars.PH,
                },
            )

    m.fs.compressor = HelmIsentropicCompressor(
            default={
                "property_package": m.fs.prop_water,
                }
            )

    return m


def main():
    m = make_model()
    return m


if __name__ == "__main__":
    m = main()
    print(degrees_of_freedom(m))
