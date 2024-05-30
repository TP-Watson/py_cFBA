__version__ = (0, 0, 0)
__all__ = []

import pandas as pd
import numpy as np
import libsbml
from optlang import Constraint, Variable, Model
from time import time


def cFBA_backbone_from_S_matrix(S_matrix):
    """
    Generate an Excel backbone for a cFBA model based on the provided Stoichiometric matrix.

    Parameters:
    S_matrix (pd.DataFrame): Stoichiometric matrix where rows represent metabolites and columns represent reactions.

    Returns:
    dict: Dictionary containing user inputs for model configuration.
    """

    # Get metabolite and reaction labels
    metabolites = list(S_matrix.index)
    reactions = list(S_matrix.columns)

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
            input(
                "Enter the indices of imbalanced metabolites (comma-separated): "
            ).split(","),
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
    use_capacities = (
        input("\nDo you want to test capacities? (yes/no): ").lower().strip()
    )
    if use_capacities == "yes":
        data["use_capacities"] = True

        # Show imbalanced metabolites for capacity testing
        print("\nWhich of the following (imbalanced metabolites) are catalysts:")
        for i, metabolite in enumerate(imbalanced_metabolites, 1):
            print(f"{i}. {metabolite}")
        catalyst_indices = list(
            map(
                int,
                input(
                    "Enter the indices of catalyst metabolites (comma-separated): "
                ).split(","),
            )
        )
        catalyst_metabolites = [imbalanced_metabolites[i - 1] for i in catalyst_indices]

        # Store catalysts in data dictionary
        data["catalysts"] = catalyst_metabolites

    return data, dt


def generate_cFBA_excel_sheet(S_matrix, data, output_file_name):
    """
    Generate an Excel backbone for a cFBA model based on the provided Stoichiometric matrix and user data.

    Parameters:
    S_matrix (pd.DataFrame): Stoichiometric matrix where rows represent metabolites and columns represent reactions.
    data (dict): Dictionary containing user inputs for model configuration.
    output_file_name (str): Name of the output Excel file.

    Returns:
    None
    """
    # Create Excel writer object
    writer = pd.ExcelWriter(output_file_name, engine="xlsxwriter")

    # Tab 1: S_mat
    pd.DataFrame(S_matrix).to_excel(writer, sheet_name="S_mat")

    # Tab 2: Imbalanced_mets
    imbalanced_mets_df = pd.DataFrame(
        {"Met": data["Imbalanced metabolites"], "w_matrix": ""}
    )
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

    # Save the Excel file
    writer.save()


def excel_to_sbml(excel_file, output_file):
    """
    Converts metabolic model data from an Excel file (format for cFBA) to SBML format and saves it to an output file.

    Parameters:
    - excel_file (str): Path to the Excel file containing metabolic model data.
    - output_file (str): Path to the output SBML file to be created.

    Returns:
    - None
    """

    # Read the Excel sheet containing stoichiometric matrix
    S_mat = pd.read_excel(excel_file, sheet_name="S_mat", header=0, index_col=0)
    S = np.array(S_mat)

    # number of metabolites and reactions
    mr, nr = S.shape

    # Reaction and metabolite labels
    rxns = list(S_mat)
    mets = list(S_mat.index)

    # Balanced and imbalanced metabolites with w matrix
    data = pd.read_excel(excel_file, sheet_name="Imbalanced_mets", header=0)
    imbalanced_mets = list(data["Met"])
    w = np.array(data["w_matrix"])

    # Read the capacity matrices from Excel
    Acap = np.array(
        pd.read_excel(excel_file, sheet_name="A_cap", header=0, index_col=None)
    )
    Bcap = np.array(
        pd.read_excel(excel_file, sheet_name="B_cap", header=0, index_col=None)
    )

    # Read time and variable bounds
    t = np.array(
        pd.read_excel(excel_file, sheet_name="lb_var", header=None, index_col=0)
    )[0]
    nt = np.size(t)
    dt = t[1] - t[0]  # Steps in time

    low_b_var = np.array(
        pd.read_excel(excel_file, sheet_name="lb_var", header=0, index_col=0)
    )
    upp_b_var = np.array(
        pd.read_excel(excel_file, sheet_name="ub_var", header=0, index_col=0)
    )

    # Get metabolite and reaction labels
    rxns = list(S_mat.columns)
    mets = list(S_mat.index)

    # Create SBML model
    document = libsbml.SBMLDocument(3, 1)
    model = document.createModel()
    model.setId("Basic_model_cFBA")
    model.setName("Basic model cFBA")

    # Define compartments
    cytoplasm = model.createCompartment()
    cytoplasm.setId("cytoplasm")
    cytoplasm.setName("Cytoplasm")
    cytoplasm.setSpatialDimensions(3)  # Optional, set spatial dimensions
    cytoplasm.setConstant(True)  # Optional, set compartment as constant

    extracellular_space = model.createCompartment()
    extracellular_space.setId("extracellular_space")
    extracellular_space.setName("Extracellular Space")
    extracellular_space.setSpatialDimensions(3)  # Optional, set spatial dimensions
    extracellular_space.setConstant(True)  # Optional

    # Define metabolites
    for metabolite_id in mets:
        metabolite = model.createSpecies()
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
        reaction = model.createReaction()
        reaction.setId(reaction_id)

        # Set required attributes
        reaction.setReversible(True)  # or False, depending on your model
        reaction.setFast(False)  # or True, depending on your model

        # Set time-specific upper and lower bounds for the reaction
        for time_index, time_point in enumerate(t):
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


def read_sbml_file(sbml_file):
    """
    Read an SBML file and return the SBML document.

    Parameters:
        sbml_file (str): Path to the SBML file.

    Returns:
        libsbml.SBMLDocument: SBML document object.
    """
    # Create an SBML reader object
    reader = libsbml.SBMLReader()

    # Read the SBML file and obtain the SBML document
    document = reader.readSBML(sbml_file)

    # Check for any errors in the SBML document
    if document.getNumErrors() > 0:
        # Print any encountered errors
        print("Encountered the following SBML errors:")
        print(document.getErrorLog().toString())

    # Return the SBML document
    return document


def parse_compartments(model):
    """
    Parse compartments from the SBML model and return compartment dictionary.

    Parameters:
        model (libsbml.Model): SBML model object.

    Returns:
        dict: Dictionary containing compartment information.
    """
    # Initialize an empty dictionary to store compartment information
    compartments = {}

    # Iterate over each compartment in the model
    for i in range(model.getNumCompartments()):
        # Get the compartment object
        compartment = model.getCompartment(i)

        # Extract compartment information and add it to the dictionary
        compartments[compartment.getId()] = {
            "size": compartment.getSize(),
            # Additional attributes can be added here if needed
        }

    # Return the dictionary containing compartment information
    return compartments


def parse_species(model):
    """
    Parse species from the SBML model and return species dictionary.

    Parameters:
        model (libsbml.Model): SBML model object.

    Returns:
        dict: Dictionary containing species information.
    """
    # Initialize an empty dictionary to store species information
    species = {}

    # Iterate over each species in the model
    for i in range(model.getNumSpecies()):
        # Get the species object
        specie = model.getSpecies(i)

        # Extract annotation and check if species is imbalanced
        annotation = specie.getAnnotationString()
        is_imbalanced = False
        w_contribution = None
        if "<isImbalanced>true</isImbalanced>" in annotation:
            is_imbalanced = True

            # Extract wContribution if available
            start_idx = annotation.find("<wContribution>")
            end_idx = annotation.find("</wContribution>")
            if start_idx != -1 and end_idx != -1:
                w_contribution = float(
                    annotation[start_idx + len("<wContribution>") : end_idx]
                )

        # Add species information to the dictionary
        species[specie.getId()] = {
            "compartment": specie.getCompartment(),
            "imbalanced": is_imbalanced,
            "w_contribution": w_contribution,
        }

    # Return the dictionary containing species information
    return species


def parse_reactions(model):
    """
    Parse reactions from the SBML model and return reaction dictionary.

    Parameters:
        model (libsbml.Model): SBML model object.

    Returns:
        dict: Dictionary containing reaction information.
    """
    # Initialize an empty dictionary to store reaction information
    reactions = {}

    # Iterate over each reaction in the model
    for i in range(model.getNumReactions()):
        # Get the reaction object
        reaction = model.getReaction(i)

        # Initialize a dictionary to store reaction data
        reaction_data = {
            "reactants": {},
            "products": {},
            "kinetic_law": {},
            "annotations": {},
        }

        # Extract reactants and their stoichiometry
        for j in range(reaction.getNumReactants()):
            reactant = reaction.getReactant(j)
            reaction_data["reactants"][
                reactant.getSpecies()
            ] = reactant.getStoichiometry()

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


def initialize_S_matrix(species, reactions):
    """
    Initialize the stoichiometry matrix S.

    Parameters:
        species (dict): Dictionary containing species data.
        reactions (dict): Dictionary containing reaction data.

    Returns:
        S (numpy.ndarray): Initialized stoichiometry matrix.
        mets (list): List of metabolite labels.
        rxns (list): List of reaction labels.
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
    for reaction_index, (reaction_id, reaction_data) in enumerate(reactions.items()):
        # Update stoichiometry for reactants
        for metabolite_id, stoichiometry in reaction_data["reactants"].items():
            metabolite_index = mets.index(metabolite_id)
            S[metabolite_index, reaction_index] -= stoichiometry
        # Update stoichiometry for products
        for metabolite_id, stoichiometry in reaction_data["products"].items():
            metabolite_index = mets.index(metabolite_id)
            S[metabolite_index, reaction_index] += stoichiometry

    return S, mets, rxns


def extract_imbalanced_metabolites(species, mets, S):
    """
    Extract indices and data for balanced and imbalanced metabolites.

    Parameters:
        species (dict): Dictionary containing species data.
        mets (list): List of metabolite labels.
        S (numpy.ndarray): Stoichiometry matrix.

    Returns:
        indices_balanced (list): Indices of balanced metabolites.
        indices_imbalanced (list): Indices of imbalanced metabolites.
        imbalanced_mets (list): List of imbalanced metabolite labels.
        balanced_mets (list): List of balanced metabolite labels.
        w (numpy.ndarray): Array of w_contribution values for imbalanced metabolites.
        Sb (numpy.ndarray): Stoichiometry matrix for balanced metabolites.
        Si (numpy.ndarray): Stoichiometry matrix for imbalanced metabolites.
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

    return (
        indices_balanced,
        indices_imbalanced,
        imbalanced_mets,
        balanced_mets,
        w,
        Sb,
        Si,
    )


def extract_kinetic_parameters(reactions):
    """
    Extract lower and upper bounds for kinetic parameters.

    Parameters:
        reactions (dict): Dictionary containing reaction data.

    Returns:
        low_b_var (numpy.ndarray): Array of lower bounds for kinetic parameters.
        upp_b_var (numpy.ndarray): Array of upper bounds for kinetic parameters.
    """
    low_b_var = []
    upp_b_var = []

    # Iterate over each reaction to extract kinetic parameters
    for reaction_id, reaction_data in reactions.items():
        kinetic_law_data = reaction_data["kinetic_law"]

        # Extract lower bounds for kinetic parameters
        lb_values = [
            kinetic_law_data.get(f"LB_{i}", 0)
            for i in range(len(kinetic_law_data))
            if f"LB_{i}" in kinetic_law_data
        ]  # Get LB_i values or default to 0
        # Extract upper bounds for kinetic parameters
        ub_values = [
            kinetic_law_data.get(f"UB_{i}", 0)
            for i in range(len(kinetic_law_data))
            if f"LB_{i}" in kinetic_law_data
        ]  # Get UB_i values or default to 0

        low_b_var.append(lb_values)
        upp_b_var.append(ub_values)

    low_b_var = np.array(low_b_var)
    upp_b_var = np.array(upp_b_var)

    return low_b_var, upp_b_var


def generate_time_components(low_b_var):
    """
    Generate time components based on the size of the lower bound array.

    Parameters:
        low_b_var (numpy.ndarray): Array of lower bounds for kinetic parameters.

    Returns:
        nt (int): Number of time steps.
    """
    # Determine the number of time steps
    nt = np.size(low_b_var[0])

    return nt


def generate_B_and_A_matrices(reactions, imbalanced_mets):
    """
    Generate B and A matrices for capacities.

    Parameters:
        reactions (dict): Dictionary containing reaction data.
        imbalanced_mets (list): List of imbalanced metabolite labels.

    Returns:
        Bcap (numpy.ndarray): B matrix.
        Acap (numpy.ndarray): A matrix.
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
            position_start = reaction_data["annotation"].find("<catalyzedBy>") + len(
                "<catalyzedBy>"
            )
            position_end = reaction_data["annotation"].find("</catalyzedBy>")
            catalyst_i = reaction_data["annotation"][position_start:position_end]

            # Find the A value
            position_start = reaction_data["annotation"].find("<A_value>") + len(
                "<A_value>"
            )
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

    return Bcap, Acap


def generate_LP_cFBA(sbml_file, quotas, dt):
    """
    Generate LP problem for constrained flux balance analysis (cFBA).

    Parameters:
        sbml_file (str): Path to the SBML file.
        quotas (list): List containing all the quota definitions for the model in the form [type, metabolite, time, value].
        dt (float): Time step increment.

    Returns:
        cons (list): List of constraints.
        Mk (numpy.ndarray): Array of metabolite amounts over time.
        imbalanced_mets (list): List of imbalanced metabolite labels.
        nm (int): Number of metabolites.
        nr (int): Number of reactions.
        nt (int): Number of time steps.
    """

    # Read SBML file and parse model components
    document = read_sbml_file(sbml_file)
    model = document.getModel()
    compartments = parse_compartments(model)
    species = parse_species(model)
    reactions = parse_reactions(model)

    # Initialize matrices and extract metabolite information
    S, mets, rxns = initialize_S_matrix(species, reactions)
    indices_balanced, indices_imbalanced, imbalanced_mets, balanced_mets, w, Sb, Si = (
        extract_imbalanced_metabolites(species, mets, S)
    )
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
            [
                Variable(
                    f"{rxns[j]}__{i}", lb=low_b_var[j, i - 1], ub=upp_b_var[j, i - 1]
                )
                for j in range(nr)
            ]
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
        cons.append(con)

    # Capacity constraints (Acap * vk <= Bcap * Mk-1)
    if Acap.shape == (0,):
        print("No catalytic capacities defined")
    else:
        for i, row in enumerate(np.dot(Acap, vk) - np.dot(Bcap, Mk[:, :-1])):
            for j, exp in enumerate(row, 1):
                con = Constraint(exp, ub=0, name=f"capacity{ i }_{ j }")
                cons.append(con)

    return cons, Mk, imbalanced_mets, nm, nr, nt


def create_lp_problem(alpha, cons_new, Mk, imbalanced_mets):
    """
    Create LP problem to optimize cyclic growth rate.

    Parameters:
        alpha (float): Initial value for cyclic growth rate.
        cons_new (list): List of constraints.
        Mk (numpy.ndarray): Array of metabolite amounts over time.
        imbalanced_mets (list): List of imbalanced metabolite labels.

    Returns:
        prob (Model): LP problem object.
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


def find_alpha(cons, Mk, imbalanced_mets):
    """
    Find the optimal value for the cyclic growth rate alpha.

    Parameters:
        cons (list): List of constraints.
        Mk (numpy.ndarray): Array of metabolite amounts over time.
        imbalanced_mets (list): List of imbalanced metabolite labels.

    Returns:
        alpha (float): Optimal value for cyclic growth rate.
        prob (Model): LP problem object after optimization.
    """
    start = time()

    alpha = 1

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

    return alpha, prob


def get_fluxes_amounts(sbml_file, prob, dt):
    """
    Obtain fluxes and metabolite amounts over time from cFBA simulations.

    Parameters:
        sbml_file (str): Path to the SBML file.
        prob (Model): LP problem object.
        dt (float): Time step increment.

    Returns:
        fluxes (numpy.ndarray): Array of fluxes over time.
        amounts (numpy.ndarray): Array of metabolite amounts over time.
        t (numpy.ndarray): Time array.
    """
    # Read SBML file and parse model components
    document = read_sbml_file(sbml_file)
    model = document.getModel()

    compartments = parse_compartments(model)
    species = parse_species(model)
    reactions = parse_reactions(model)

    # Initialize matrices and extract relevant data
    S, mets, rxns = initialize_S_matrix(species, reactions)
    indices_balanced, indices_imbalanced, imbalanced_mets, balanced_mets, w, Sb, Si = (
        extract_imbalanced_metabolites(species, mets, S)
    )
    low_b_var, upp_b_var = extract_kinetic_parameters(reactions)
    nt = generate_time_components(low_b_var)
    t = np.arange(0, nt * dt, dt)

    # Initialize arrays to store fluxes and amounts
    fluxes = np.zeros((len(rxns), nt - 1))
    amounts = np.zeros((len(imbalanced_mets), nt))

    # Extract fluxes and amounts from LP problem variables
    for var in prob.variables.values():
        name, j = var.name.split("__")
        j = int(j) - 1
        try:
            i = rxns.index(name)
            fluxes[i, j] = var.primal
        except ValueError:
            i = imbalanced_mets.index(name)
            amounts[i, j + 1] = var.primal

    # Calculate metabolite amounts based on stoichiometry and fluxes
    amounts[:, 1:] = np.dot(Si, fluxes) * dt
    amounts = np.cumsum(amounts, axis=1)

    return fluxes, amounts, t
