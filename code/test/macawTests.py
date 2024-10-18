import cobra
from macaw.main import dead_end_test, duplicate_test

model = cobra.io.load_yaml_model("model/Human-GEM.yml")
(dead_end_results, dead_end_edges) = dead_end_test(model)
(duplicate_results, duplicate_edges) = duplicate_test(model)
output = dead_end_results.merge(duplicate_results)
output.to_csv('data/testResults/macaw_results.csv', index = False)