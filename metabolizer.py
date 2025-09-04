# %%
#Import packages
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd

# %%
# User inputs
Ngen = 3  # Number of generations
reaction_library = pd.DataFrame({
    'Scheme_name': ['Hydrogenolysis', 'Vicinal Dehalogenation'],
    'Reaction_expression': ['[#6;A:1][#17,#35,#53]>>[#6;A:1]', '[#17,#35,#53][#6;A:1][#6;A:2][#17,#35,#53]>>[#6;A:1]=[#6;A:2]'],
    'Rank': [4, 4],
    'Reactivity_rule': [True, True],  # Simplified for demonstration
    'Selectivity_rule': [False, False]  # Simplified for demonstration
})

# %%
# Initialize parent molecule
mol = ['C(Cl)CBr']  # Parent SMILES
gen = [0]
ParentID = [0]
F_deg = [0.0]
production = [1.0]
accumulation = [0.0]
N_prod = [1]
mol_index = 0

# Initialize reaction objects
reactions = [AllChem.ReactionFromSmarts(exp) for exp in reaction_library['Reaction_expression']]

# %%
# Main loop for generations
for i in range(1, Ngen + 1):
    Np = 0
    # Loop over products in previous generation as parents in current generation
    for p in range(N_prod[i - 1]):
        j = mol_index + p
        # Ensure F_deg is long enough
        if j >= len(F_deg):
            F_deg.append(0.0)
        else:
            F_deg[j] = 0.0
        # Assign the parent ID for bookkeeping
        parent_idx = j
        reactant_mol = Chem.MolFromSmiles(mol[j])
        for k in range(len(reactions)):
            reaction = reactions[k]
            if reaction_library['Reactivity_rule'][k]:  # Simplified rule check
                F = 7 ** reaction_library['Rank'][k]
                num_exp_products = reaction.GetNumProductTemplates()
                # Reactions in Metabolizer libraries have only one reactant; define single element tuple for this reactant
                reactant_tuple = (reactant_mol,)
                try:
                    products = reaction.RunReactants(reactant_tuple)
                except Exception:
                    continue
                num_prod = len(products)
                for product_set in products:
                    F_deg[parent_idx] += F
                    for product in product_set:
                        product_smiles = Chem.MolToSmiles(product)
                        Np += 1
                        # Extend lists for new molecule
                        mol.append(product_smiles)
                        gen.append(i)
                        ParentID.append(parent_idx)
                        production.append(F)   # Store formation for production calculation
                        accumulation.append(0.0)  # Placeholder for accumulation calculation
    N_prod.append(Np)
    mol_index += N_prod[i - 1]
    prod_index = mol_index + N_prod[i]

# %%
# Calculate production values
for j in range(1, prod_index):
    production[j] = production[j]*production[ParentID[j]] / F_deg[ParentID[j]]

# Calculate accumulation values
for j in range(1, prod_index):
    accumulation[j] = production[j] - production[ParentID[j]] * F_deg[j] / F_deg[ParentID[j]]

# Output results
for idx, smiles in enumerate(mol):
    print(f"Molecule {idx}: SMILES={smiles}, Generation={gen[idx]}, Production={production[idx]}, Accumulation={accumulation[idx]}")


