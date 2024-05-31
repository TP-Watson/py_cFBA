"""Core functionality of the Python cFBA Toolbox."""

__all__ = []

from time import time
from typing import Any, cast

import libsbml
import numpy as np
import optlang.interface
import pandas as pd
from numpy.typing import NDArray

from py_cfba.typing import (
    AlphaOutput,
    CapacityMatrices,
    Constraint,
    ExtractedImbalancedMetabolites,
    FileName,
    FluxOutput,
    InitSMatrix,
    KineticParamsBounds,
    LPProblemOutput,
    Model,
    Quota,
    ReactionData,
    ReactionDict,
    SpeciesData,
    SpeciesDict,
    Variable,
)


def cFBA_backbone_from_S_matrix(S_matrix: pd.DataFrame) -> tuple[dict[str, Any], float]:
    """
    Generate an Excel backbone for a cFBA model based on the provided Stoichiometric matrix.

    Args:
        S_matrix: Stoichiometric matrix where rows represent metabolites and columns represent reactions.

    Returns:
        dict: Dictionary containing user inputs for model configuration.
        dt: Time gap
    """

    # Get metabolite and reaction labels
    metabolites = list(S_matrix.index)

    # Initialize data dictionary to store user inputs
    data = {
        "Imbalanced metabolites": [],
        "total_time": None,
        "dt": None,
        "use_capacities": False,
        "catalysts": {},
    }

    # Ask user to select imbalanced metabolites
    print("\n-------------- Imbalanced metabolites --------------")
    print("Select the imbalanced metabolites:")
    for i, metabolite in enumerate(metabolites, 1):
        print(f"{i}. {metabolite}")
    imbalanced_indices = list(
        map(
            int,
            input("Enter the indices of imbalanced metabolites (comma-separated): ").split(","),
        )
    )
    imbalanced_metabolites = [metabolites[i - 1] for i in imbalanced_indices]

    # Store imbalanced metabolites in data dictionary
    data["Imbalanced metabolites"] = imbalanced_metabolites

    # Ask for total time and dt
    print("\n-------------- Simulation time --------------")
    data["total_time"] = float(input("\nEnter the total time for simulation: "))
    data["dt"] = float(input("Enter the time gap (dt) to be simulated: "))
    dt = data["dt"]

    # Ask if user wants to test capacities
    print("\n-------------- Enzyme capacities --------------")
    use_capacities = input("\nDo you want to test capacities? (yes/no): ").lower().strip()
    if use_capacities == "yes":
        data["use_capacities"] = True

        # Show imbalanced metabolites for capacity testing
        print("\nWhich of the following (imbalanced metabolites) are catalysts:")
        for i, metabolite in enumerate(imbalanced_metabolites, 1):
            print(f"{i}. {metabolite}")
        catalyst_indices = list(
            map(
                int,
                input("Enter the indices of catalyst metabolites (comma-separated): ").split(","),
            )
        )
        catalyst_metabolites = [imbalanced_metabolites[i - 1] for i in catalyst_indices]

        # Store catalysts in data dictionary
        data["catalysts"] = catalyst_metabolites

    return data, dt


def generate_cFBA_excel_sheet(S_matrix: pd.DataFrame, data: dict[str, Any], output_file_name: FileName) -> None:
    """
    Generate an Excel backbone for a cFBA model based on the provided Stoichiometric matrix and user data.

    Args:
        S_matrix: Stoichiometric matrix where rows represent metabolites and columns represent reactions.
        data: Dictionary containing user inputs for model configuration.
        output_file_name: Name of the output Excel file.
    """
    # Create Excel writer object
    with pd.ExcelWriter(output_file_name, engine="xlsxwriter") as writer:

        # Tab 1: S_mat
        pd.DataFrame(S_matrix).to_excel(writer, sheet_name="S_mat")

        # Tab 2: Imbalanced_mets
        imbalanced_mets_df = pd.DataFrame({"Met": data["Imbalanced metabolites"], "w_matrix": ""})
        imbalanced_mets_df.to_excel(writer, sheet_name="Imbalanced_mets", index=False)

        # Tab 3: lb_var
        # Calculate time points
        time_points = np.arange(0, data["total_time"] + data["dt"], data["dt"])
        decimal_places = len(str(data["dt"]).split(".")[1])
        # Round each value to the determined number of decimal places
        time_points_rounded = [round(value, decimal_places) for value in time_points]
        lb_var_df = pd.DataFrame(0, index=S_matrix.columns, columns=time_points_rounded)
        lb_var_df.to_excel(writer, sheet_name="lb_var")

        # Tab 4: ub_var
        ub_var_df = pd.DataFrame(1000, index=S_matrix.columns, columns=time_points_rounded)
        ub_var_df.to_excel(writer, sheet_name="ub_var")

        # Tab 5: A_cap

        a_cap_df = pd.DataFrame(index=S_matrix.columns).T
        a_cap_df.to_excel(writer, sheet_name="A_cap", index=False)

        # Tab 6: B_cap
        if data["use_capacities"]:
            # Determine the number of catalyzers
            num_catalyzers = len(data["catalysts"])

            # Create DataFrame filled with zeros
            b_cap_df = pd.DataFrame(
                0,
                index=data["Imbalanced metabolites"],
                columns=list(range(1, num_catalyzers + 1)),
            )

            # Fill 1s in the corresponding positions
            for i, catalyzer in enumerate(data["catalysts"], 1):
                b_cap_df.loc[catalyzer, i] = 1
            b_cap_df.T.to_excel(writer, sheet_name="B_cap", index=False)
        else:
            b_cap_df = pd.DataFrame()
            b_cap_df.to_excel(writer, sheet_name="B_cap", index=False)


def excel_to_sbml(excel_file: FileName, output_file: FileName) -> None:
    """
    Converts metabolic model data from an Excel file (format for cFBA) to SBML format and saves it to an output file.

    Args:
        excel_file: Path to the Excel file containing metabolic model data.
        output_file: Path to the output SBML file to be created.
    """

    # Read the Excel sheet containing stoichiometric matrix
    S_mat = pd.read_excel(excel_file, sheet_name="S_mat", header=0, index_col=0)
    S = np.array(S_mat)

    # Reaction and metabolite labels
    rxns = list(S_mat)
    mets = list(S_mat.index)

    # Balanced and imbalanced metabolites with w matrix
    data = pd.read_excel(excel_file, sheet_name="Imbalanced_mets", header=0)
    imbalanced_mets = list(data["Met"])
    w = np.array(data["w_matrix"])

    # Read the capacity matrices from Excel
    Acap = np.array(pd.read_excel(excel_file, sheet_name="A_cap", header=0, index_col=None))
    Bcap = np.array(pd.read_excel(excel_file, sheet_name="B_cap", header=0, index_col=None))

    # Read time and variable bounds
    t = np.array(pd.read_excel(excel_file, sheet_name="lb_var", header=None, index_col=0))[0]

    low_b_var = np.array(pd.read_excel(excel_file, sheet_name="lb_var", header=0, index_col=0))
    upp_b_var = np.array(pd.read_excel(excel_file, sheet_name="ub_var", header=0, index_col=0))

    # Get metabolite and reaction labels
    rxns = list(S_mat.columns)
    mets = list(S_mat.index)

    # Create SBML model
    document = libsbml.SBMLDocument(3, 1)
    model: libsbml.Model = document.createModel()
    model.setId("Basic_model_cFBA")
    model.setName("Basic model cFBA")

    # Define compartments
    cytoplasm: libsbml.Compartment = model.createCompartment()
    cytoplasm.setId("cytoplasm")
    cytoplasm.setName("Cytoplasm")
    cytoplasm.setSpatialDimensions(3)  # Optional, set spatial dimensions
    cytoplasm.setConstant(True)  # Optional, set compartment as constant

    extracellular_space: libsbml.Compartment = model.createCompartment()
    extracellular_space.setId("extracellular_space")
    extracellular_space.setName("Extracellular Space")
    extracellular_space.setSpatialDimensions(3)  # Optional, set spatial dimensions
    extracellular_space.setConstant(True)  # Optional

    # Define metabolites
    for metabolite_id in mets:
        metabolite: libsbml.Species = model.createSpecies()
        metabolite.setId(metabolite_id)
        metabolite.setCompartment("cytoplasm")  # Set compartment

        # Set required attributes
        metabolite.setBoundaryCondition(False)
        metabolite.setHasOnlySubstanceUnits(True)
        metabolite.setConstant(False)  # can change dynamnically

        # Add annotation to indicate if metabolite is imbalanced
        if metabolite_id in imbalanced_mets:
            annotation = f"<annotation><metadata><isImbalanced>true</isImbalanced><wContribution>{w[imbalanced_mets.index(metabolite_id)]}</wContribution></metadata></annotation>"
        else:
            annotation = "<annotation><metadata><isImbalanced>false</isImbalanced></metadata></annotation>"

        metabolite.setAnnotation(annotation)

    # Define reactions
    for reaction_id, stoichiometry in zip(rxns, S.T):
        reaction: libsbml.Reaction = model.createReaction()
        reaction.setId(reaction_id)

        # Set required attributes
        reaction.setReversible(True)  # or False, depending on your model
        reaction.setFast(False)  # or True, depending on your model

        # Set time-specific upper and lower bounds for the reaction
        for time_index, _ in enumerate(t):
            lower_bound = float(low_b_var[rxns.index(reaction_id), time_index])
            upper_bound = float(upp_b_var[rxns.index(reaction_id), time_index])

            # Create kinetic law if it doesn't exist
            if not reaction.isSetKineticLaw():
                reaction.createKineticLaw()

            # Create parameters for lower and upper bounds
            lb_parameter = reaction.getKineticLaw().createParameter()
            lb_parameter.setId(f"LB_{time_index}")
            lb_parameter.setValue(lower_bound)
            lb_parameter.setConstant(True)

            ub_parameter = reaction.getKineticLaw().createParameter()
            ub_parameter.setId(f"UB_{time_index}")
            ub_parameter.setValue(upper_bound)
            ub_parameter.setConstant(True)

        # Check if reaction is catalyzed by imbalanced metabolites based on A_cap matrix
        A_column = Acap[:, rxns.index(reaction_id)]
        for i, A_value in enumerate(A_column):
            if A_value != 0:
                # Find the imbalanced metabolite catalyzing this reaction based on B_cap matrix
                B_row = Bcap[i]
                imb_met_index = np.where(B_row == 1)[0][0]
                imb_met_id = imbalanced_mets[imb_met_index]

                # Get the 1/kcat value from A_cap matrix
                inv_kcat_value = A_value

                # Create annotation for catalysis by imbalanced metabolite
                annotation = f"<annotation><metadata><catalyzedBy>{imb_met_id}</catalyzedBy><A_value>{inv_kcat_value}</A_value></metadata></annotation>"
                reaction.setAnnotation(annotation)

        # Add reactants and products
        for met_id, stoich_coeff in zip(mets, stoichiometry):
            if stoich_coeff != 0:
                if stoich_coeff < 0:
                    reactant = reaction.createReactant()
                    reactant.setSpecies(met_id)
                    reactant.setStoichiometry(float(abs(stoich_coeff)))
                    reactant.setConstant(True)
                else:
                    product = reaction.createProduct()
                    product.setSpecies(met_id)
                    product.setStoichiometry(float(stoich_coeff))
                    product.setConstant(True)

    # Write SBML document to file
    libsbml.writeSBMLToFile(document, output_file)
    print(
        f"SBML document with metabolites information and catalysis annotations has been created and saved to {output_file}."
    )


def read_sbml_file(sbml_file: FileName) -> libsbml.SBMLDocument:
    """
    Read an SBML file and return the SBML document.

    Args:
        sbml_file: Path to the SBML file.

    Returns:
        document: SBML document object.
    """
    sbml_file_ = cast(Any, sbml_file)

    # Create an SBML reader object
    reader = libsbml.SBMLReader()

    # Read the SBML file and obtain the SBML document
    document: libsbml.SBMLDocument = reader.readSBML(sbml_file_)

    # Check for any errors in the SBML document
    if document.getNumErrors() > 0:
        # Print any encountered errors
        print("Encountered the following SBML errors:")
        print(document.getErrorLog().toString())

    # Return the SBML document
    return document


def parse_compartments(sbml_model: libsbml.Model) -> dict[str, dict[str, Any]]:
    """
    Parse compartments from the SBML model and return compartment dictionary.

    Args:
        model: SBML model object.

    Returns:
        compartments: Dictionary containing compartment information.
    """
    # Initialize an empty dictionary to store compartment information
    compartments = {}

    # Iterate over each compartment in the model
    for i in range(sbml_model.getNumCompartments()):
        # Get the compartment object
        compartment = sbml_model.getCompartment(i)

        # Extract compartment information and add it to the dictionary
        compartments[compartment.getId()] = {
            "size": compartment.getSize(),
            # Additional attributes can be added here if needed
        }

    # Return the dictionary containing compartment information
    return compartments


def parse_species(sbml_model: libsbml.Model) -> SpeciesDict:
    """
    Parse species from the SBML model and return species dictionary.

    Args:
        model: SBML model object.

    Returns:
        species: Dictionary containing species information.
    """
    # Initialize an empty dictionary to store species information
    species: SpeciesDict = {}

    # Iterate over each species in the model
    for i in range(sbml_model.getNumSpecies()):
        # Get the species object
        metabolite: libsbml.Species = sbml_model.getSpecies(i)

        # Extract annotation and check if species is imbalanced
        annotation = metabolite.getAnnotationString()
        is_imbalanced = False
        w_contribution = None
        if "<isImbalanced>true</isImbalanced>" in annotation:
            is_imbalanced = True

            # Extract wContribution if available
            start_idx = annotation.find("<wContribution>")
            end_idx = annotation.find("</wContribution>")
            if start_idx != -1 and end_idx != -1:
                w_contribution = float(annotation[start_idx + len("<wContribution>") : end_idx])

        # Add species information to the dictionary
        species_data: SpeciesData = {
            "compartment": metabolite.getCompartment(),
            "imbalanced": is_imbalanced,
            "w_contribution": w_contribution,
        }
        species[metabolite.getId()] = species_data

    # Return the dictionary containing species information
    return species


def parse_reactions(sbml_model: libsbml.Model) -> ReactionDict:
    """
    Parse reactions from the SBML model and return reaction dictionary.

    Args:
        model: SBML model object.

    Returns:
        reactions: Dictionary containing reaction information.
    """
    # Initialize an empty dictionary to store reaction information
    reactions: ReactionDict = {}

    # Iterate over each reaction in the model
    for i in range(sbml_model.getNumReactions()):
        # Get the reaction object
        reaction: libsbml.Reaction = sbml_model.getReaction(i)

        # Initialize a dictionary to store reaction data
        reaction_data: ReactionData = {
            "reactants": {},
            "products": {},
            "kinetic_law": {},
            "annotation": "",
        }

        # Extract reactants and their stoichiometry
        for j in range(reaction.getNumReactants()):
            reactant = reaction.getReactant(j)
            reaction_data["reactants"][reactant.getSpecies()] = reactant.getStoichiometry()

        # Extract products and their stoichiometry
        for j in range(reaction.getNumProducts()):
            product = reaction.getProduct(j)
            reaction_data["products"][product.getSpecies()] = product.getStoichiometry()

        # Extract kinetic law parameters if available
        kinetic_law = reaction.getKineticLaw()
        if kinetic_law:
            for j in range(kinetic_law.getNumParameters()):
                parameter = kinetic_law.getParameter(j)
                reaction_data["kinetic_law"][parameter.getId()] = parameter.getValue()

        # Get annotation for the reaction
        annotation = reaction.getAnnotationString()
        reaction_data["annotation"] = annotation

        # Add reaction data to the reactions dictionary
        reactions[reaction.getId()] = reaction_data

    # Return the dictionary containing reaction information
    return reactions


def initialize_S_matrix(species: SpeciesDict, reactions: ReactionDict) -> InitSMatrix:
    """
    Initialize the stoichiometry matrix S.

    Args:
        species: Dictionary containing species data.
        reactions: Dictionary containing reaction data.

    Returns:
        S: Initialized stoichiometry matrix.
        mets: List of metabolite labels.
        rxns: List of reaction labels.
    """
    # Extract metabolite and reaction labels
    mets = list(species.keys())
    rxns = list(reactions.keys())

    # Get the number of metabolites and reactions
    nm = len(mets)
    nr = len(rxns)

    # Initialize the stoichiometry matrix S with zeros
    S = np.zeros((nm, nr))

    # Populate the S matrix based on reactants and products of each reaction
    for reaction_index, (_, reaction_data) in enumerate(reactions.items()):
        # Update stoichiometry for reactants
        for metabolite_id, stoichiometry in reaction_data["reactants"].items():
            metabolite_index = mets.index(metabolite_id)
            S[metabolite_index, reaction_index] -= stoichiometry
        # Update stoichiometry for products
        for metabolite_id, stoichiometry in reaction_data["products"].items():
            metabolite_index = mets.index(metabolite_id)
            S[metabolite_index, reaction_index] += stoichiometry

    return InitSMatrix(S, mets, rxns)


def extract_imbalanced_metabolites(species: SpeciesDict, mets: list[str], S: NDArray) -> ExtractedImbalancedMetabolites:
    """
    Extract indices and data for balanced and imbalanced metabolites.

    Args:
        species: Dictionary containing species data.
        mets: List of metabolite labels.
        S: Stoichiometry matrix.

    Returns:
        indices_balanced: Indices of balanced metabolites.
        indices_imbalanced: Indices of imbalanced metabolites.
        imbalanced_mets: List of imbalanced metabolite labels.
        balanced_mets: List of balanced metabolite labels.
        w: Array of w_contribution values for imbalanced metabolites.
        Sb: Stoichiometry matrix for balanced metabolites.
        Si: Stoichiometry matrix for imbalanced metabolites.
    """
    indices_balanced = []
    indices_imbalanced = []
    imbalanced_mets = []
    w = []

    # Iterate over species data to categorize metabolites
    for metabolite_id, metabolite_data in species.items():
        if metabolite_data["imbalanced"]:
            indices_imbalanced.append(mets.index(metabolite_id))
            imbalanced_mets.append(metabolite_id)
            w_contribution = species[metabolite_id]["w_contribution"]
            w.append(w_contribution)
        else:
            indices_balanced.append(mets.index(metabolite_id))

    # Create a list of balanced metabolites
    balanced_mets = [met for met in mets if met not in imbalanced_mets]
    w = np.array(w)

    # Extract stoichiometric matrices for balanced and imbalanced metabolites
    Sb = S[indices_balanced, :]
    Si = S[indices_imbalanced, :]

    return ExtractedImbalancedMetabolites(
        indices_balanced,
        indices_imbalanced,
        imbalanced_mets,
        balanced_mets,
        w,
        Sb,
        Si,
    )


def extract_kinetic_parameters(reactions: ReactionDict) -> KineticParamsBounds:
    """
    Extract lower and upper bounds for kinetic parameters.

    Args:
        reactions: Dictionary containing reaction data.

    Returns:
        low_b_var: Array of lower bounds for kinetic parameters.
        upp_b_var: Array of upper bounds for kinetic parameters.
    """
    low_b_var = []
    upp_b_var = []

    # Iterate over each reaction to extract kinetic parameters
    for _, reaction_data in reactions.items():
        kinetic_law_data = reaction_data["kinetic_law"]

        # Extract lower bounds for kinetic parameters
        lb_values = [
            kinetic_law_data.get(f"LB_{i}", 0) for i in range(len(kinetic_law_data)) if f"LB_{i}" in kinetic_law_data
        ]  # Get LB_i values or default to 0
        # Extract upper bounds for kinetic parameters
        ub_values = [
            kinetic_law_data.get(f"UB_{i}", 0) for i in range(len(kinetic_law_data)) if f"LB_{i}" in kinetic_law_data
        ]  # Get UB_i values or default to 0

        low_b_var.append(lb_values)
        upp_b_var.append(ub_values)

    low_b_var = np.array(low_b_var)
    upp_b_var = np.array(upp_b_var)

    return KineticParamsBounds(low_b_var, upp_b_var)


def generate_time_components(low_b_var: NDArray) -> int:
    """
    Generate time components based on the size of the lower bound array.

    Args:
        low_b_var: Array of lower bounds for kinetic parameters.

    Returns:
        nt: Number of time steps.
    """
    # Determine the number of time steps
    return np.size(low_b_var, axis=1)


def generate_B_and_A_matrices(reactions: ReactionDict, imbalanced_mets: list[str]) -> CapacityMatrices:
    """
    Generate B and A matrices for capacities.

    Args:
        reactions: Dictionary containing reaction data.
        imbalanced_mets: List of imbalanced metabolite labels.

    Returns:
        Bcap: B matrix.
        Acap: A matrix.
    """
    # Extract reaction IDs
    rxns = list(reactions.keys())

    # Initialize dictionaries to store catalysts and A values
    catalyzed_by = {}
    A_values = {}

    # Iterate over each reaction in the reactions dictionary
    for reaction_id, reaction_data in reactions.items():
        # Check if the reaction contains catalyst information
        if reaction_data["annotation"].find("<catalyzedBy>") > 0:
            # Find the catalyst
            position_start = reaction_data["annotation"].find("<catalyzedBy>") + len("<catalyzedBy>")
            position_end = reaction_data["annotation"].find("</catalyzedBy>")
            catalyst_i = reaction_data["annotation"][position_start:position_end]

            # Find the A value
            position_start = reaction_data["annotation"].find("<A_value>") + len("<A_value>")
            position_end = reaction_data["annotation"].find("</A_value>")
            A_val_i = reaction_data["annotation"][position_start:position_end]

            # Store the catalyst and A value in the dictionaries
            catalyzed_by[reaction_id] = imbalanced_mets.index(catalyst_i)
            A_values[reaction_id] = A_val_i

    # Generate the B matrix (imb_mets x unique_catalyzers)
    unique_catalyzer = list(set(catalyzed_by.values()))
    Bcap = np.eye(len(imbalanced_mets))[unique_catalyzer, :]

    # Generate the A matrix (unique_catalyzers x reactions)
    Acap = np.zeros([Bcap.shape[0], len(rxns)])
    for i, arr in enumerate(Bcap):
        catalyst_i = imbalanced_mets[np.where(arr == 1)[0][0]]
        for reaction_id, cat_pos in catalyzed_by.items():
            reaction_pos = rxns.index(reaction_id)
            if catalyst_i == imbalanced_mets[cat_pos]:
                Acap[i, reaction_pos] = A_values[reaction_id]

    return CapacityMatrices(Bcap, Acap)


def generate_LP_cFBA(sbml_file: FileName, quotas: list[Quota], dt: float) -> LPProblemOutput:
    """
    Generate LP problem for constrained flux balance analysis (cFBA).

    Args:
        sbml_file: Path to the SBML file.
        quotas: List containing all the quota definitions for the model in the form [type, metabolite, time, value].
        dt: Time step increment.

    Returns:
        cons: List of constraints.
        Mk: Array of metabolite amounts over time.
        imbalanced_mets: List of imbalanced metabolite labels.
        nm: Number of metabolites.
        nr: Number of reactions.
        nt: Number of time steps.
    """

    # Read SBML file and parse model components
    document = read_sbml_file(sbml_file)
    model = document.getModel()
    species = parse_species(model)
    reactions = parse_reactions(model)

    # Initialize matrices and extract metabolite information
    S, mets, rxns = initialize_S_matrix(species, reactions)
    _, _, imbalanced_mets, balanced_mets, w, Sb, Si = extract_imbalanced_metabolites(species, mets, S)
    low_b_var, upp_b_var = extract_kinetic_parameters(reactions)
    nt = generate_time_components(low_b_var)
    Bcap, Acap = generate_B_and_A_matrices(reactions, imbalanced_mets)

    nm = len(mets)
    nr = len(rxns)

    # Linear Programming (LP) problem setup
    cons = []

    # Define variables: fluxes and starting amounts
    vk = np.array(
        [
            [Variable(f"{rxns[j]}__{i}", lb=low_b_var[j, i - 1], ub=upp_b_var[j, i - 1]) for j in range(nr)]
            for i in range(1, nt)
        ]
    ).T
    M0 = np.array([[Variable(f"{c}__0", lb=0, ub=1000) for c in imbalanced_mets]]).T

    # Calculate metabolite amounts over time
    Mk = np.dot(Si, vk) * dt
    Mk = np.hstack((M0, Mk))
    Mk = np.cumsum(Mk, axis=1)

    # Non-negative amounts constraint (Mk >= 0)
    for j, c in enumerate(imbalanced_mets):
        for i in range(1, nt):
            con = Constraint(Mk[j, i], lb=0, ub=1000, name=f"{c}__{i}")
            cons.append(con)

    # Steady-state mass balances constraint (Sb*vk = 0)
    for i, row in enumerate(np.dot(Sb, vk)):
        for j, exp in enumerate(row, 1):
            con = Constraint(exp, lb=0, ub=0, name=f"{balanced_mets[i]}_{j}")
            cons.append(con)

    # Starting biomass constraint (wT*M0 = 1)
    con = Constraint(np.dot(w.T, Mk[:, 0]), lb=1, ub=1, name="starting_biomass")
    cons.append(con)

    # Quota constraints
    for i, entry in enumerate(quotas):
        # Extract each quota info
        cons_type, met, time_point, value = entry
        # B_mat: position of metabolite
        B = np.zeros((1, len(imbalanced_mets)))
        B[0, imbalanced_mets.index(met)] = 1
        # Lelft hand side equation
        exp = np.dot(B, Mk[:, time_point]) - value
        # Create constraint based on its type
        if cons_type == "equality":
            con = Constraint(exp[0], lb=0, ub=0, name=f"quota_eq{i}")
        elif cons_type == "min":
            con = Constraint(exp[0], lb=0, name=f"quota_min{i}")
        elif cons_type == "max":
            con = Constraint(exp[0], ub=0, name=f"quota_max{i}")
        else:
            raise ValueError("Unrecognized type of constraint. Valid types are 'equality', 'min', and 'max'.")
        cons.append(con)

    # Capacity constraints (Acap * vk <= Bcap * Mk-1)
    if Acap.shape == (0,):
        print("No catalytic capacities defined")
    else:
        for i, row in enumerate(np.dot(Acap, vk) - np.dot(Bcap, Mk[:, :-1])):
            for j, exp in enumerate(row, 1):
                con = Constraint(exp, ub=0, name=f"capacity{ i }_{ j }")
                cons.append(con)

    return LPProblemOutput(cons, Mk, imbalanced_mets, nm, nr, nt)


def create_lp_problem(
    alpha: float, cons_new: list[optlang.interface.Constraint], Mk: NDArray, imbalanced_mets: list[str]
) -> optlang.interface.Model:
    """
    Create LP problem to optimize cyclic growth rate.

    Args:
        alpha: Initial value for cyclic growth rate.
        cons_new: List of constraints.
        Mk: Array of metabolite amounts over time.
        imbalanced_mets: List of imbalanced metabolite labels.

    Returns:
        prob: LP problem object.
    """
    # Initialize LP problem
    prob = Model()

    # Add cyclic growth constraints
    for i, exp in enumerate(Mk[:, -1] - Mk[:, 0] * alpha):
        con = Constraint(exp, lb=0, ub=0, name=f"alpha_{imbalanced_mets[i]}")
        cons_new.append(con)

    # Add constraints to the LP problem
    prob.add(cons_new)

    return prob


def find_alpha(cons: list[optlang.interface.Constraint], Mk: NDArray, imbalanced_mets: list[str]) -> AlphaOutput:
    """
    Find the optimal value for the cyclic growth rate alpha.

    Args:
        cons: List of constraints.
        Mk: Array of metabolite amounts over time.
        imbalanced_mets: List of imbalanced metabolite labels.

    Returns:
        alpha: Optimal value for cyclic growth rate.
        prob: LP problem object after optimization.
    """
    start = time()

    alpha = 1.0

    # Iterate to find the upper bound for alpha
    while True:
        prob = create_lp_problem(2 * alpha, [*cons], Mk, imbalanced_mets)
        status = prob.optimize()
        if status == "optimal":
            alpha *= 2
        else:
            break

    elapsed_time = time() - start

    start = time()

    delta = alpha / 2

    # Use binary search to find the optimal alpha
    while delta > 1e-10:
        prob = create_lp_problem(alpha + delta, [*cons], Mk, imbalanced_mets)
        status = prob.optimize()
        if status == "optimal":
            alpha += delta
        delta = delta / 2

    prob = create_lp_problem(alpha - delta, [*cons], Mk, imbalanced_mets)
    prob.optimize()

    elapsed_time = (time() - start) / 60
    print(f"{elapsed_time:.2f} min")

    return AlphaOutput(alpha, prob)


def get_fluxes_amounts(sbml_file: FileName, prob: optlang.interface.Model, dt: float) -> FluxOutput:
    """
    Obtain fluxes and metabolite amounts over time from cFBA simulations.

    Args:
        sbml_file: Path to the SBML file.
        prob: LP problem object.
        dt: Time step increment.

    Returns:
        fluxes: Array of fluxes over time.
        amounts: Array of metabolite amounts over time.
        t: Time array.
    """
    # Read SBML file and parse model components
    document = read_sbml_file(sbml_file)
    model: libsbml.Model = document.getModel()

    species = parse_species(model)
    reactions = parse_reactions(model)

    # Initialize matrices and extract relevant data
    S, mets, rxns = initialize_S_matrix(species, reactions)
    ext = extract_imbalanced_metabolites(species, mets, S)
    low_b_var, _ = extract_kinetic_parameters(reactions)
    nt = generate_time_components(low_b_var)
    t = np.arange(0, nt * dt, dt)

    # Initialize arrays to store fluxes and amounts
    fluxes = np.zeros((len(rxns), nt - 1))
    amounts = np.zeros((len(ext.imbalanced_mets), nt))

    # Extract fluxes and amounts from LP problem variables
    for var in prob.variables.values():
        name, j = var.name.split("__")
        j = int(j) - 1
        try:
            i = rxns.index(name)
            fluxes[i, j] = var.primal
        except ValueError:
            i = ext.imbalanced_mets.index(name)
            amounts[i, j + 1] = var.primal

    # Calculate metabolite amounts based on stoichiometry and fluxes
    amounts[:, 1:] = np.dot(ext.Si, fluxes) * dt
    amounts = np.cumsum(amounts, axis=1)

    return FluxOutput(fluxes, amounts, t)
