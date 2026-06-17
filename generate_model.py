from auxin_signalling_models.AuxinSignallingPathway import AuxinSignallingPathway


def main():
    model = AuxinSignallingPathway(nARFs=1, nIAAs=1, include_arf_transcription=True)

    model.add_state("ATHB8", category="vasular")
    model.add_state("X", category="vascular")


    model.add_state("D__ATHB8__ATHB8", category="dimer")
    model.add_state("D__ATHB8__X", category="dimer")
    model.add_state("D__X__X", category="dimer")

    ## mRNA state variables
    model.add_state("R__ATHB8", category="mRNA")
    model.add_state("R__X", category="mRNA")

    ## Add protein production/decay rates
    model.add_model_parameter("k__ATHB8", 1.0, "production rate for ATHB8",
                              "production"
                              )

    model.add_model_parameter("k__X", 1.0, "production rate for X",
                              "production"
                              )

    model.add_model_parameter("d__ATHB8", 1.0, "decay rate for ATHB8",
                              "decay"
                              )

    model.add_model_parameter("d__X", 1.0, "decay rate for X",
                              "decay"
                              )

    model.add_reaction([model.empty_set_node.name], ["ATHB8"],
                       ["k__ATHB8"], ["d__ATHB8"])

    model.add_reaction([model.empty_set_node.name], ["X"],
                       ["k__X"], ["d__X"])

    ## Dimerisation reactions
    model.add_model_parameter("k__D__ATHB8__ATHB8", 1.0,
                              "hetero-dimerisation rate for ATHB8", "dimerisation")
    model.add_model_parameter("d__D__ATHB8__ATHB8", 1.0,
                              "decay rate for ATHB8 heterodimers", "dimerisation")

    model.add_reaction(["ATHB8", "ATHB8"], ["D__ATHB8__ATHB8"],
                       fwd_rates=["k__D__ATHB8__ATHB8"],
                       bwd_rates=["d__D__ATHB8__ATHB8"])

    model.add_model_parameter("k__D__ATHB8__X", 1.0,
                              "Dimerisation rate for ATHB8-X", "dimerisation")
    model.add_model_parameter("d__D__ATHB8__X", 1.0,
                              "decay rate for ATHB8-X dimers", "dimerisation")


    model.add_reaction(["ATHB8", "X"], ["D__ATHB8__X"],
                       fwd_rates=["k__D__ATHB8__X"],
                       bwd_rates=["d__D__ATHB8__X"])

    model.add_model_parameter("k__D__X__X", 1.0,
                              "hetero-dimerisation rate for X", "dimerisation")
    model.add_model_parameter("d__D__X__X", 1.0,
                              "decay rate for X heterodimers", "dimerisation")

    model.add_reaction(["X", "X"], ["D__X__X"],
                       fwd_rates=["k__D__X__X"],
                       bwd_rates=["d__D__X__X"])


    ## Promoter binding reactions (using a separate promoter)
    model.add_state("G2", category="promoter")

    model.add_state("G2__D__ATHB8__ATHB8", category="promoter")
    model.add_state("G2__D__ATHB8__X", category="promoter")
    model.add_state("G2__D__X__X", category="promoter")
    model.add_state("G2__D__ATHB8__ATHB8", category="promoter")

    model.add_model_parameter("k__G2__D__ATHB8__ATHB8", category="promoter binding")
    model.add_model_parameter("d__G2__D__ATHB8__ATHB8", category="promoter binding")
    model.add_reaction(["G2"], ["G2__D__ATHB8__ATHB8"], ["k__G2__D__ATHB8__ATHB8"], [],
                       enzymes=["D__ATHB8__ATHB8"])
    model.add_reaction(["G2"], ["G2__D__ATHB8__ATHB8"], [], ["d__G2__D__ATHB8__ATHB8"],
                       enzymes=["D__ATHB8__ATHB8"])

    model.add_model_parameter("k__G2__D__ATHB8__X", category="promoter binding")
    model.add_model_parameter("d__G2__D__ATHB8__X", category="promoter binding")
    model.add_reaction(["G2"], ["G2__D__ATHB8__X"], ["k__G2__D__ATHB8__X"], [])
    model.add_reaction(["G2"], ["G2__D__ATHB8__X"], [], ["d__G2__D__ATHB8__X"])

    model.add_model_parameter("k__G2__D__X__X", category="promoter binding")
    model.add_model_parameter("d__G2__D__X__X", category="promoter binding")
    model.add_reaction(["G2"], ["G2__D__X__X"], ["k__G2__D__X__X"], [])
    model.add_reaction(["G2"], ["G2__D__X__X"], [], ["d__G2__D__X__X"])

    ## TODO add transcription rates
    for mrna in ["R__ATHB8", "R__X"]:
        for promoter_state in ["G2", "G2__D__ATHB8__ATHB8", "G2__D_ATHB8_X_", "G2__D__X__X"]:
            t_rate = model.add_model_parameter(f"k__{promoter_state}__{mrna}", 1.0,
                                               f"Transcription rate for {mrna} when G2 promoter in state {promoter_state}",
                                               category="transcription")
            model.add_reaction([promoter_state], [mrna], [t_rate])

    print(model.generate_pretty_print_state_dict())
    print(model.generate_pretty_print_parameter_dict())


if __name__ == "__main__":
    main()
