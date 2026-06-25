import os

import numpy as np
import sympy as sp
from sympy.printing.julia import julia_code

from typing import Callable, List
from types import MethodType

from auxin_signalling_models.AuxinSignallingPathway import AuxinSignallingPathway
from auxin_signalling_models.ReducedAuxinSignallingPathway\
    import ReducedAuxinSignallingPathway

import matplotlib.pyplot as plt


def main():
    model = generate_model()

    output_dir = os.path.join("output")
    os.makedirs(output_dir, exist_ok=True)
    model.draw_graph()
    plt.savefig(os.path.join(output_dir, "test_graph"))
    plt.clf()


def generate_model():
    model = ATHB8Model(nARFs=0, nIAAs=0,
                       include_arf_transcription=True)

    return model


class ATHB8Model(ReducedAuxinSignallingPathway):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.add_state("ATHB8", category="vasular")
        self.add_state("X", category="vascular")

        self.add_state("D__ATHB8__ATHB8", category="dimer")
        self.add_state("D__ATHB8__X", category="dimer")
        self.add_state("D__X__X", category="dimer")

        ## mRNA state variables
        self.add_state("R__ATHB8", category="mRNA")
        self.add_state("R__X", category="mRNA")

        ## Add protein production/decay rates
        self.add_model_parameter("k__ATHB8", 1.0e-2, "production rate for ATHB8",
                                "production"
                                )

        self.add_model_parameter("k__X", 1.0e-2, "production rate for X",
                                "production"
                                )

        self.add_model_parameter("d__ATHB8", 1.0, "decay rate for ATHB8",
                                "decay"
                                )

        self.add_model_parameter("d__X", 1.0, "decay rate for X",
                                "decay"
                                )

        self.add_reaction([self.empty_set_node.name], ["ATHB8"],
                          [], ["d__ATHB8"])

        self.add_reaction([self.empty_set_node.name], ["ATHB8"],
                          ["k__ATHB8"], [], enzymes=["R__ATHB8"])

        self.add_reaction([self.empty_set_node.name], ["X"],
                          [], ["d__X"])

        self.add_reaction([self.empty_set_node.name], ["X"],
                          ["k__X"], [], enzymes=["R__X"])

        ## Dimerisation reactions
        self.add_model_parameter("k__D__ATHB8__ATHB8", 1.0e-1,
                                 "hetero-dimerisation rate for ATHB8", "dimerisation")
        self.add_model_parameter("d__D__ATHB8__ATHB8", 1.0e-1,
                                 "decay rate for ATHB8 heterodimers", "dimerisation")

        self.add_reaction(["ATHB8", "ATHB8"], ["D__ATHB8__ATHB8"],
                          fwd_rates=["k__D__ATHB8__ATHB8"],
                          bwd_rates=["d__D__ATHB8__ATHB8"])

        self.add_model_parameter("k__D__ATHB8__X", 1.0e-1,
                                 "Dimerisation rate for ATHB8-X", "dimerisation")
        self.add_model_parameter("d__D__ATHB8__X", 1.0e-1,
                                 "decay rate for ATHB8-X dimers", "dimerisation")

        self.add_reaction(["ATHB8", "X"], ["D__ATHB8__X"],
                          fwd_rates=["k__D__ATHB8__X"],
                          bwd_rates=["d__D__ATHB8__X"])

        self.add_model_parameter("k__D__X__X", 1.0,
                                 "hetero-dimerisation rate for X", "dimerisation")
        self.add_model_parameter("d__D__X__X", 1.0,
                                 "decay rate for X heterodimers", "dimerisation")

        self.add_reaction(["X", "X"], ["D__X__X"],
                          fwd_rates=["k__D__X__X"],
                          bwd_rates=["d__D__X__X"])

        ## Promoter binding reactions (using a separate promoter)
        self.add_state("G2", category="promoter")

        self.add_state("G2__D__ATHB8__ATHB8", category="promoter")
        self.add_state("G2__D__ATHB8__X", category="promoter")
        self.add_state("G2__D__X__X", category="promoter")

        self.add_model_parameter("k__G2__D__ATHB8__ATHB8",
                                 category="promoter binding")
        self.add_model_parameter("d__G2__D__ATHB8__ATHB8",
                                 category="promoter binding")
        self.add_reaction(["G2"], ["G2__D__ATHB8__ATHB8"],
                          ["k__G2__D__ATHB8__ATHB8"], [],
                          enzymes=["D__ATHB8__ATHB8"])

        self.add_reaction(["G2"], ["G2__D__ATHB8__ATHB8"], [],
                          ["d__G2__D__ATHB8__ATHB8"])

        self.add_model_parameter("k__G2__D__ATHB8__X", category="promoter binding")
        self.add_model_parameter("d__G2__D__ATHB8__X", category="promoter binding")
        self.add_reaction(["G2"], ["G2__D__ATHB8__X"], ["k__G2__D__ATHB8__X"], [],
                          enzymes=["D__ATHB8__X"])
        self.add_reaction(["G2"], ["G2__D__ATHB8__X"], [], ["d__G2__D__ATHB8__X"],
                          )

        self.add_model_parameter("k__G2__D__X__X", category="promoter binding")
        self.add_model_parameter("d__G2__D__X__X", category="promoter binding")
        self.add_reaction(["G2"], ["G2__D__X__X"], ["k__G2__D__X__X"], [],
                          enzymes=["D__X__X"])
        self.add_reaction(["G2"], ["G2__D__X__X"], [], ["d__G2__D__X__X"])

        for mrna in ["R__ATHB8", "R__X"]:

            self.add_reaction([mrna], [self.empty_set_node.name], ["d__R"])

            for promoter_state in ["G2", "G2__D__ATHB8__ATHB8",
                                   "G2__D__ATHB8__X", "G2__D__X__X"]:
                t_rate = self.add_model_parameter(f"k__{promoter_state}__{mrna}",
                                                  0.0,
                                                  f"Transcription rate for {mrna}"
                                                  "when G2 promoter in state"
                                                  f"{promoter_state}",
                                                  category="transcription")
                self.add_reaction([self.empty_set_node.name], [mrna],
                                  [t_rate.name])

        print(self.generate_pretty_print_state_dict())
        print(self.generate_pretty_print_parameter_dict())

        self.add_model_parameter("k__G__ARF_1__ARF_1__R__ATHB8", 1.0, "",
                                 "transcription")
        self.add_model_parameter("k__G__ARF_1__ARF_1__R__X", 1.0, "",
                                 "transcription")

        self.add_reaction([self.empty_set_node.name], ["R__ATHB8"],
                          ["k__G__ARF_1__ARF_1__R__ATHB8"],
                          enzymes=["G__ARF_1__ARF_1"])
        self.add_reaction([self.empty_set_node.name], ["R__X"],
                          ["k__G__ARF_1__ARF_1__R__X"],
                          enzymes=["G__ARF_1__ARF_1"])


if __name__ == "__main__":
    main()
