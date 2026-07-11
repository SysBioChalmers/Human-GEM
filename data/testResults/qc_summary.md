#### MACAW: dead-end and duplicate tests

```
Starting dead-end test...
 - Found 1380 dead-end metabolites.
 - Found 1137 reactions incapable of sustaining steady-state fluxes in either direction due to these dead-ends.
 - Found 1369 reversible reactions that can only carry steady-state fluxes in a single direction due to dead-ends.
Starting duplicate test...
 - Skipping redox duplicates because no redox_pairs and/or proton_ids were provided.
 - Found 377 reactions that were some type of duplicate:
   - 0 were completely identical to at least one other reaction.
   - 13 involve the same metabolites but go in the opposite direction or have the opposite reversibility as at least one other reaction.
   - 377 involve the same metabolites but with different coefficients as at least one other reaction.
```

#### Mass and charge balance

```
Unbalanced reactions (excluding boundary and biomass): 277 (87 mass, 240 charge)
```
