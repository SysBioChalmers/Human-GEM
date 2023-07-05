import cobra
import memote.support.basic

if __name__ == "__main__":
    try:
        model = cobra.io.load_yaml_model("model/Human-GEM.yml")

        # Use the MEMOTE way of testing
        # https://github.com/opencobra/memote/blob/c8bd3fe75e2d955deaec824d1e93ebe60e754710/tests/test_for_support/test_for_basic.py#L908-L909
        rxns, num = memote.support.basic.find_duplicate_reactions(model)
        assert num == 0, "Duplicate reactions found: {}".format(rxns)

        # Make sure all reactions have at least one metabolite
        for r in model.reactions:
            assert len(r.metabolites) > 0, "Reaction with no metabolite found: {}".format(r.id)
    except AssertionError as e:
        print(e)
        exit(1)