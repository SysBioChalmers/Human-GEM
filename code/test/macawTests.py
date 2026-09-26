import cobra
from macaw.main import dead_end_test, duplicate_test

model = cobra.io.load_yaml_model("model/Human-GEM.yml")
(dead_end_results, dead_end_edges) = dead_end_test(model)
(duplicate_results, duplicate_edges) = duplicate_test(model)
output = dead_end_results.merge(duplicate_results)

# Only reactions with a finding are written. "only when going ..." marks a reversible
# reaction limited to one direction, which can still carry flux, so it is not a finding.
duplicate_columns = [c for c in output.columns if c.startswith("duplicate_test")]
dead_end = ~output["dead_end_test"].eq("ok") & ~output["dead_end_test"].str.startswith("only when going")
duplicate = ~output[duplicate_columns].isin(["ok", "N/A"]).all(axis=1)
output[dead_end | duplicate].to_csv("data/testResults/macaw_results.tsv", sep="\t", index=False)
