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

    # Why "forward" solver?
    ic_dict = model.get_default_initial_conditions()
    ic_dict = {k: 1.0 for k in ic_dict}
    ts = np.linspace(0, 3000, 2)
    param_dict = model.get_default_parameters()
    param_dict = {k: 1.0 for k in param_dict}
    solver = model.make_forward_solver()

    res = solver(ts, param_dict, ic_dict, 0.99)

    print(res)


def generate_model():
    model = ATHB8Model(nARFs=1, nIAAs=0,
                       include_arf_transcription=False)

    return model


class ATHB8Model(ReducedAuxinSignallingPathway):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.add_state("ATHB8", category="vasular")
        self.add_state("X", category="vascular")

        self.add_state("D__ATHB8__ATHB8", category="dimer")
        # self.add_state("D__ATHB8__X", category="dimer")
        self.add_state("D__X__X", category="dimer")

        # mRNA state variables
        self.add_state("R__ATHB8", category="mRNA")
        self.add_state("R__X", category="mRNA")

        # Add protein production/decay rates
        self.add_model_parameter("k__ATHB8", 1.0e-1, "production rate for ATHB8",
                                 "production"
                                 )

        self.add_model_parameter("k__X", 1.0e-1, "production rate for X",
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

        # RNA decay
        self.add_reaction(["R__X"], [self.empty_set_node.name], ["d__R"])
        self.add_reaction(["R__ATHB8"], [self.empty_set_node.name], ["d__R"])

        # Dimerisation reactions
        self.add_model_parameter("k__D__ATHB8__ATHB8", 1.0e-1,
                                 "hetero-dimerisation rate for ATHB8",
                                 "dimerisation")
        self.add_model_parameter("d__D__ATHB8__ATHB8", 1.0e-1,
                                 "decay rate for ATHB8 heterodimers",
                                 "dimerisation")

        self.add_reaction(["ATHB8", "ATHB8"], ["D__ATHB8__ATHB8"],
                          fwd_rates=["k__D__ATHB8__ATHB8"],
                          bwd_rates=["d__D__ATHB8__ATHB8"])

        self.add_model_parameter("k__D__ATHB8__X", 1.0e-1,
                                 "Dimerisation rate for ATHB8-X", "dimerisation")
        self.add_model_parameter("d__D__ATHB8__X", 1.0e-1,
                                 "decay rate for ATHB8-X dimers", "dimerisation")

        # self.add_reaction(["ATHB8", "X"], ["D__ATHB8__X"],
        #                   fwd_rates=["k__D__ATHB8__X"],
        #                   bwd_rates=["d__D__ATHB8__X"])

        self.add_model_parameter("k__D__X__X", 1.0,
                                 "hetero-dimerisation rate for X", "dimerisation")
        self.add_model_parameter("d__D__X__X", 1.0,
                                 "decay rate for X heterodimers", "dimerisation")

        self.add_reaction(["X", "X"], ["D__X__X"],
                          fwd_rates=["k__D__X__X"],
                          bwd_rates=["d__D__X__X"])

        ## Promoter binding reactions (using a separate promoter)
        self.add_state("G2", category="promoter")

        self.add_state("G2__ATHB8__ATHB8", category="promoter")
        # self.add_state("G2__ATHB8__X", category="promoter")
        self.add_state("G2__X__X", category="promoter")

        self.add_model_parameter("k__G2__ATHB8__ATHB8",
                                 category="promoter binding")
        self.add_model_parameter("d__G2__ATHB8__ATHB8",
                                 category="promoter binding")
        self.add_reaction(["G2"], ["G2__ATHB8__ATHB8"],
                          ["k__G2__ATHB8__ATHB8"], [],
                          enzymes=["D__ATHB8__ATHB8"])

        self.add_reaction(["G2"], ["G2__ATHB8__ATHB8"], [],
                          ["d__G2__ATHB8__ATHB8"])

        # self.add_model_parameter("k__G2__ATHB8__X", category="promoter binding")
        # self.add_model_parameter("d__G2__ATHB8__X", category="promoter binding")
        # self.add_reaction(["G2"], ["G2__ATHB8__X"], ["k__G2__ATHB8__X"], [],
        #                   enzymes=["D__ATHB8__X"])
        # self.add_reaction(["G2"], ["G2__ATHB8__X"], [], ["d__G2__ATHB8__X"],
        #                   )

        self.add_model_parameter("k__G2__X__X", category="promoter binding")
        self.add_model_parameter("d__G2__X__X", category="promoter binding")
        self.add_reaction(["G2"], ["G2__X__X"], ["k__G2__X__X"], [],
                          enzymes=["D__X__X"])
        self.add_reaction(["G2"], ["G2__X__X"], [], ["d__G2__X__X"])

        for mrna in ["R__ATHB8", "R__X"]:

            self.add_reaction([mrna], [self.empty_set_node.name], ["d__R"])

            promoter_states =  ["G2",
                                "G2__ATHB8__ATHB8",
                                # "G2__ATHB8__X",
                                "G2__X__X"]
 
            for promoter_state in promoter_states:
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

        self.add_state("G2__ARF_1__ARF_1", category="promoter")
        self.add_model_parameter("k__G2__ARF_1__ARF_1",
                                 1.0, "Binding rate for G2 promoter with ARF-ARF dimer")

        self.add_model_parameter("d__G2__ARF_1__ARF_1",
                                 1.0, "Unbindnig rate for G2 promoter with ARF-ARF dimer")

        self.add_reaction(["G2"], ["G2__ARF_1__ARF_1"],
                          ["k__G2__ARF_1__ARF_1"],
                          [], enzymes=["D__ARF_1__ARF_1"])

        self.add_reaction(["G2__ARF_1__ARF_1"],
                          ["G2"],
                          ["d__G2__ARF_1__ARF_1"],
                          [])

        self.add_model_parameter("k__G2__ARF_1__ARF_1__R__ATHB8", 1.0, "",
                                 "transcription")
        self.add_model_parameter("k__G2__ARF_1__ARF_1__R__X", 1.0, "",
                                 "transcription")

        self.add_reaction([self.empty_set_node.name], ["R__ATHB8"],
                          ["k__G2__ARF_1__ARF_1__R__ATHB8"],
                          enzymes=["G2__ARF_1__ARF_1"])
        self.add_reaction([self.empty_set_node.name], ["R__X"],
                          ["k__G2__ARF_1__ARF_1__R__X"],
                          enzymes=["G2__ARF_1__ARF_1"])

        # self.add_state("G2__ARF_1__iaa_1", category="promoter")
        # self.add_model_parameter("k__G2__ARF_1__iaa_1",
        #                          0.1, "Binding rate for G2 promoter with ARF-Aux/IAA dimer")

        # self.add_model_parameter("d__G2__ARF_1__iaa_1",
        #                          0.1, "Unbindnig rate for G2 promoter with ARF-Aux/IAA dimer")

        # self.add_reaction(["G2"], ["G2__ARF_1__iaa_1"],
        #                   ["k__G2__ARF_1__iaa_1"],
        #                   [], enzymes=["D__ARF_1__iaa_1"])

        # self.add_reaction(["G2__ARF_1__iaa_1"],
        #                   ["G2"],
        #                   ["d__G2__ARF_1__iaa_1"],
        #                   [], [])

        # self.add_model_parameter("k__G2__ARF_1__iaa_1__R__ATHB8", 1.0, "",
        #                          "transcription")
        # self.add_model_parameter("k__G2__ARF_1__iaa_1__R__X", 1.0, "",
        #                          "transcription")

        # self.add_reaction([self.empty_set_node.name], ["R__ATHB8"],
        #                   ["k__G2__ARF_1__iaa_1__R__ATHB8"],
        #                   enzymes=["G2__ARF_1__iaa_1"])
        # self.add_reaction([self.empty_set_node.name], ["R__X"],
        #                   ["k__G2__ARF_1__iaa_1__R__X"],
        #                   enzymes=["G2__ARF_1__iaa_1"])

        states_to_remove = ["G", "G__ARF_1__ARF_1"]
        for state in states_to_remove:
            self.remove_state_variable_by_name(state)

        print(list(self.graph.nodes()))


if __name__ == "__main__":
    main()
