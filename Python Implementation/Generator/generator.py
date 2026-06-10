##########################################
#           Generator                   #
#      Juliana Hirt, 2026-06-09         #    
##########################################


##########################################
#         IMPORT STATEMENTS.             #
########################################## 

import random # used for random sampling and choices
import uuid # used for generating unique identifiers for each person in the synthetic population
import json # used for loading and parsing JSON schema files
import pandas as pd # used for data manipulation and creating the synthetic population dataframe
import numpy as np # used for numerical operations and random number generation
import yaml # used for loading and parsing YAML configuration files
import warnings # used to suppress warnings during data generation and manipulation
warnings.filterwarnings("ignore")

##########################################
#         GLOBAL VARIABLES.              #
##########################################
config = None # will hold the configuration settings loaded from the YAML file
schema = None # will hold the JSON schema loaded from the specified file
schema_columns = None # will hold the flattened schema metadata as a list of dictionaries
schema_df = None # will hold the flattened schema metadata as a pandas DataFrame
data_df = None # will hold the generated synthetic population data as a pandas DataFrame
rng = None # will hold the random number generator instance initialized with the seed from the configuration


##########################################
#              FUNCTIONS                 #
########################################## 

def read_config(config_path: str) -> dict:
    """
    Reads a YAML configuration file and returns its contents as a dictionary.
    """

    with open(config_path, "r", encoding="utf-8") as file:
        return yaml.safe_load(file)
    

def load_json_schema(path: str) -> dict:
    """
    Loads a JSON schema from the specified file path and returns it as a dictionary.
    """
    with open(path, "r", encoding="utf-8") as file:
        return json.load(file)
    
def json_type_to_pandas_dtype(json_type):
    """
    Convert a JSON Schema type to a Pandas dtype.
    """

    if isinstance(json_type, list):
        json_type = [t for t in json_type if t != "null"]
        json_type = json_type[0] if json_type else "object"

    type_mapping = {
        "string": "string",
        "integer": "Int64",
        "number": "float64",
        "boolean": "boolean",
        "array": "object",
        "object": "object",
    }

    return type_mapping.get(json_type, "object")


def flatten_json_schema(schema: dict, parent_key: str = "") -> list[dict]:
    """
    Recursively flatten a nested JSON Schema into table-style metadata.
    """

    rows = []

    properties = schema.get("properties", {})
    required_fields = set(schema.get("required", []))

    for field_name, field_schema in properties.items():
        full_name = f"{parent_key}.{field_name}" if parent_key else field_name

        field_type = field_schema.get("type")

        # If this is a nested object, recurse into its properties
        if field_type == "object" and "properties" in field_schema:
            rows.extend(flatten_json_schema(field_schema, full_name))

        else:
            rows.append(
                {
                    "column": full_name,
                    "json_type": field_type,
                    "pandas_dtype": json_type_to_pandas_dtype(field_type),
                    "required": field_name in required_fields,
                    "description": field_schema.get("description"),
                    "enum": field_schema.get("enum"),
                    "minimum": field_schema.get("minimum"),
                    "format": field_schema.get("format"),
                }
            )

    return rows

def generate_guid() -> str:
    """
    Generates a unique identifier (GUID) using the uuid4 method.
    """
    return str(uuid.uuid4())


def generate_age_str(age_config: dict) -> str:
    """
    Generates an age category based on the age configuration.
    """
    return rng.choice(
        age_config["categories"],
        p=age_config["weights"]
    )


def generate_comorbidity(comorbidity_config: dict) -> bool:
    """
    Generates a boolean value indicating the presence of comorbidities based on the comorbidity configuration.
    """

    return bool(
        rng.choice(
            [True, False],
            p=[
                comorbidity_config["probabilityTrue"],
                1 - comorbidity_config["probabilityTrue"]
            ]
        )
    )


def generate_social_activity(social_activity_config: dict) -> str:
    """
    Generates a social activity level based on the social activity configuration.
    """
    return rng.choice(
        social_activity_config["categories"],
        p=social_activity_config["weights"]
    )


def generate_geography(geography_config: dict) -> str:
    """
    Generates a geographic location category based on the geography configuration.
    """
    return rng.choice(
        geography_config["categories"],
        p=geography_config["weights"]
    )


def generate_mobility(mobility_config: dict) -> str:
    """ 
    Generates a mobility level category based on the mobility configuration.
    """
    return rng.choice(
        mobility_config["categories"],
        p=mobility_config["weights"]
    )


def generate_vaccine_acceptance(vaccine_acceptance_config: dict) -> bool:
    """
    Generates a boolean value indicating the presence of vaccine acceptance based on the vaccine acceptance configuration.
    """
    return bool(
        rng.choice(
            [True, False],
            p=[
                vaccine_acceptance_config["probabilityTrue"],
                1 - vaccine_acceptance_config["probabilityTrue"]
            ]
        )
    )


def get_risk_multiplier(field_name: str, field_value) -> float:
    """
    Retrieves the infection risk multiplier for a given field name and value from the configuration.
    """
    if isinstance(field_value, bool):
        lookup_value = str(field_value).lower()
    else:
        lookup_value = field_value

    return infection_risk_multipliers[field_name][lookup_value]

def set_rng_seed(seed: int) -> None:
    """
    Sets the seed for the random number generator to ensure reproducibility.
    """
    global rng
    rng = np.random.default_rng(seed)
    return rng

def infect_initial_population():
    """
    Infects a specified number of individuals in the synthetic population at the start of the simulation.
    """
    global data_df, config

    numb_to_infect = config["initialInfected"]

    infection_indices = rng.choice(data_df.index, size=numb_to_infect, replace=False).tolist()

    for idx in infection_indices:
        data_df.at[idx, "dynamic.sirvStatus"] = "I"


def update_positions(x, y, L, sig2, rng):
    """
    Update x and y positions using Brownian motion
    with periodic boundary conditions.
    """

    N = len(x)
    sig = np.sqrt(sig2)

    dx = rng.normal(loc=0, scale=sig, size=N)
    dy = rng.normal(loc=0, scale=sig, size=N)

    x = np.mod(x + dx, L)
    y = np.mod(y + dy, L)

    return x, y


def jiggle_positions(df, rho, sig2, rng):
    """
    Apply one Brownian motion jiggle to the dataframe positions.

    This function is meant to be imported and reused.
    It does not depend on global config or global data_df.
    """

    df = df.copy()

    x = df["dynamic.currentLocation.xcor"].to_numpy()
    y = df["dynamic.currentLocation.ycor"].to_numpy()

    N = len(df)
    L = np.sqrt(N / rho)

    x, y = update_positions(
        x=x,
        y=y,
        L=L,
        sig2=sig2,
        rng=rng
    )

    df["dynamic.currentLocation.xcor"] = x
    df["dynamic.currentLocation.ycor"] = y

    return df
        

def setup():
    """
    Sets up the generator by loading the configuration and schema, and preparing the dataframe structure.
    """
    global config, schema, schema_columns, schema_df, data_df, rng, infection_risk_multipliers

    # Load configuration and schema
    config = read_config("../Config/config_sample.yml")
    schema = load_json_schema("../Generator/person_template.json")

    # Flatten the JSON schema into a list of column metadata
    schema_columns = flatten_json_schema(schema)

    # Create a DataFrame from the flattened schema metadata
    schema_df = pd.DataFrame(schema_columns)

    # Set up a the dataframe structure based on the schema
    columns = schema_df["column"].tolist()
    dtypes = schema_df.set_index("column")["pandas_dtype"].to_dict()
    data_df = pd.DataFrame(columns=columns).astype(dtypes)

    #Set up risk multipliers from the configuration
    infection_risk_multipliers = config["infectionRiskMultipliers"]

    # Initialize random number generator with the seed from the configuration
    set_rng_seed(config["seed"])

def generate_synthetic_population():
    """
    Generates synthetic population data based on the configuration and schema.
    """
    global data_df, config

    populationSize = config["populationSize"]
    age_config = config["syntheticPopulation"]["static"]["age"]
    comorbidity_config = config["syntheticPopulation"]["static"]["comorbidity"]
    social_activity_config = config["syntheticPopulation"]["static"]["socialActivity"]
    geography_config = config["syntheticPopulation"]["static"]["geography"]
    mobility_config = config["syntheticPopulation"]["static"]["mobility"]
    vaccine_acceptance_config = config["syntheticPopulation"]["static"]["vaccineAcceptance"]
    infection_risk_multipliers = config["infectionRiskMultipliers"]

    rows = []

    for _ in range(populationSize):
        age = generate_age_str(age_config)
        comorbidity = generate_comorbidity(comorbidity_config)
        social_activity = generate_social_activity(social_activity_config)
        geography = generate_geography(geography_config)
        mobility = generate_mobility(mobility_config)
        vaccine_acceptance = generate_vaccine_acceptance(vaccine_acceptance_config)


        row = {
            "static.guid": generate_guid(),

            "static.age": age,
            "static.ageRiskMultiplier": get_risk_multiplier("age", age),

            "static.comorbidity": comorbidity,
            "static.comorbidityRiskMultiplier": get_risk_multiplier("comorbidity", comorbidity),

            "static.socialActivity": social_activity,
            "static.socialActivityRiskMultiplier": get_risk_multiplier("socialActivity", social_activity),

            "static.geography": geography,
            "static.geographyRiskMultiplier": get_risk_multiplier("geography", geography),

            "static.mobility": mobility,
            "static.mobilityRiskMultiplier": get_risk_multiplier("mobility", mobility),

            "static.vaccineAcceptance": vaccine_acceptance,
            "static.vaccineAcceptanceRiskMultiplier": get_risk_multiplier("vaccineAcceptance", vaccine_acceptance),

            "dynamic.vaccineStatus": False,
            "dynamic.proactiveVaccine": False,
            "dynamic.numberOfInfections": 0,
            "dynamic.sirvStatus": "S",
            "dynamic.infectedDays": 0,
            "dynamic.vaccinatedDays": 0,
            "dynamic.recoveredDays": 0,
            "dynamic.currentLocation.xcor": 0.0,
            "dynamic.currentLocation.ycor": 0.0
        }

        rows.append(row)

    print(f"Generated {len(rows)} synthetic population records.")
    print(f"Population size: {populationSize}")

    data_df = pd.DataFrame(rows)

    infect_initial_population()
    print(f"Infected initial population: {config['initialInfected']} individuals.")

    data_df.to_csv("synthetic_population2.csv", index=False)


def intialization() -> None:
    """
    Main function to execute the synthetic population generation process.
    """

    #initial set up to load config, schema, and prepare dataframe
    setup()
    print("Configuration and schema loaded successfully.")

    #generate the synthetic population data
    generate_synthetic_population()
    print("Synthetic population generated successfully.")

    #set initial positions with a random jiggle to simulate movement
    jiggle_positions(data_df, config["rho"], config["sig2"], rng)
    print("Initial positions jiggle applied successfully.")

    return data_df
    
    


if __name__ == "__main__":
    intialization()
    

