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


def generate_model(nARFs=0, nIAAs=0, **kwargs):
    model = ATHB8Model(nARFs=nARFs, nIAAs=nIAAs,
                       **kwargs)

    return model


class ATHB8Model(ReducedAuxinSignallingPathway):

    def __init__(self, nARFs=0, nIAAs=0, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.add_state("ATHB8", category="vasular")
        self.add_state("X", category="vascular")

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

        # Promoter binding reactions (using a separate promoter)
        self.add_state("G2", category="promoter")

        if nARFs == 0:

            self.add_model_parameter("_ARF_1_", 1.0, "ARF promoter binding rate", category="promoter")
            self.add_model_parameter("d_ARF_1_", 1.0, "ARF promoter unbinding rate",
                                    category="promoter")

            self.add_state("G2__ARF_1", category="promoter")
            self.add_reaction(["G2"], ["G2__ARF_1"], ["_ARF_1_"], ["d_ARF_1_"])

        else:
            for state in self.dimer_variables:
                dimer_name = state.name
                promoter_state_name = "G2__" + dimer_name[3:]

                if "iaa_1__iaa_1" in promoter_state_name:
                    continue

                state = self.add_state(promoter_state_name,
                                       category="promoter",
                                       components=[dimer_name])

                bind_var = self.add_model_parameter(
                    f"k__{promoter_state_name}",
                    1.0,
                    f"promoter binding rate for {dimer_name} and G2",
                )
                unbind_var = self.add_model_parameter(
                    f"d__{promoter_state_name}",
                    1.0,
                    f"promoter unbinding rate for {dimer_name} and G2",
                )
                self.add_reaction(["G2"], [promoter_state_name],
                                  [], [unbind_var.name]
                                  )

                self.add_reaction(["G2"], [promoter_state_name],
                                  [bind_var.name], [],
                                  enzymes=[dimer_name]
                                  )

                target_genes = ["X", "ATHB8"]
                for target_gene in target_genes:
                    transcription_rate_var_name = f"k__{promoter_state_name}__R__{target_gene}"
                    self.add_model_parameter(transcription_rate_var_name,
                                             1.0,
                                             "transcription rate")
                    self.add_reaction([], [f"R__{target_gene}"],
                                      [transcription_rate_var_name],
                                      [],
                                      enzymes=[promoter_state_name]
                                      )

        self.add_state("G2__ATHB8", category="promoter")
        # self.add_state("G2__ATHB8__X", category="promoter")
        self.add_state("G2__X", category="promoter")

        self.add_model_parameter("k__G2__ATHB8",
                                 category="promoter binding")
        self.add_model_parameter("d__G2__ATHB8",
                                 category="promoter binding")
        self.add_reaction(["G2"], ["G2__ATHB8"],
                          ["k__G2__ATHB8"], [],
                          enzymes=["ATHB8"])

        self.add_reaction(["G2"], ["G2__ATHB8"], [],
                          ["d__G2__ATHB8"])

        self.add_model_parameter("k__G2__X", 1.0, category="promoter binding")
        self.add_model_parameter("d__G2__X", 1.0, category="promoter binding")
        self.add_reaction(["G2"], ["G2__X"], ["k__G2__X"], [],
                          enzymes=["X"])
        self.add_reaction(["G2"], ["G2__X"], [], ["d__G2__X"])

        for mrna in ["R__ATHB8", "R__X"]:

            self.add_reaction([mrna], [self.empty_set_node.name], ["d__R"])

            promoter_states =  ["G2",
                                "G2__ATHB8",
                                "G2__X"]

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

        if nARFs == 0:
            self.add_model_parameter("k__G2__ARF_1",
                                     1.0, "Binding rate for G2 promoter with ARF")

            self.add_model_parameter("d__G2__ARF_1",
                                     1.0, "Unbindnig rate for G2 promoter with ARF")

            self.add_model_parameter("k__G2__ARF_1__R__ATHB8", 1.0, "",
                                     "transcription")
            self.add_model_parameter("k__G2__ARF_1__R__X", 1.0, "",
                                     "transcription")

            self.add_reaction([self.empty_set_node.name], ["R__ATHB8"],
                              ["k__G2__ARF_1__R__ATHB8"],
                              enzymes=["G2__ARF_1"])
            self.add_reaction([self.empty_set_node.name], ["R__X"],
                              ["k__G2__ARF_1__R__X"],
                              enzymes=["G2__ARF_1"])

        states_to_remove = []

        if nARFs == 0 and nIAAs == 0:
            states_to_remove = ["G"]
        for state in states_to_remove:
            self.remove_state_variable_by_name(state)


        # Add extra promoter state for 2nd binding site
        self.add_state("G2__ARF_1__ARF_1_2", category="promoter")

        bind_rate = self.add_model_parameter("k__G2__ARF_1__ARF_1_2", 1.0,
                                               category="promoter")

        unbind_rate = self.add_model_parameter("d__G2__ARF_1__ARF_1_2", 1.0,
                                               category="promoter")
        transcription_rate = self.add_model_parameter("k__G2__ARF_1__ARF_1_2__R__ATHB8",
                                                      1.0, category="transcription")

        self.add_reaction(["G2__ARF_1__ARF_1"], ["G2__ARF_1__ARF_1_2"], [bind_rate.name],
                          [], enzymes=["D__ARF_1__ARF_1"])

        self.add_reaction(["G2__ARF_1__ARF_1"], ["G2__ARF_1__ARF_1_2"], [],
                          [unbind_rate.name])

        self.add_reaction([self.empty_set_node.name], ["R__ATHB8"],
                          [transcription_rate.name], [],
                          enzymes=["G2__ARF_1__ARF_1_2"])

if __name__ == "__main__":
    main()
