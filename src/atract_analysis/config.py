from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
RAW_SOURCE_SHEET = "Feuil1"
RAW_FALLBACK_SHEETS = [
    "Queries in empty cells in red",
    "Feuil1",
]

RAW_COLUMNS = {
    "date_exam": "Date examen",
    "physician": "Physician",
    "operator_simple": "Operateur simplifie",
    "age": "Âge",
    "sex": "Gender (1 male, 2 woman)",
    "asa": "ASA",
    "anticoagulants": "Anticoaggulants",
    "antiplatelets": "Antiaggrégants",
    "location_simple": "Localisation simplifiee",
    "location_raw": "Localisation",
    "major_diameter": "Major diameter (mm)",
    "minor_diameter": "Minor Diameter (mm)",
    "surface": "Surface",
    "lesion_type": "Type de lésion (1 Polype, 2 LST Granulaire, 3 LST non granulaire)",
    "macronodule": "Macronodule",
    "jnet": "JNET (1= JNET  I, 2= JNET IIa, 3= JNET IIB, 4= JNET III)",
    "conecct": "CONECCT (1=1H, 2=Is, 3=IIA, 4=IIC, 5=IIC+, 6=III, 7=OE)",
    "traction": "TRACTION (0= non, 1 = oui)",
    "atract": "ATRACT (0= non, 1= oui)",
    "duration": "Duration of procedure Procédure (min)",
    "speed": "Vitesse",
    "fibrosis": "Fibrosis",
    "perforation": "Intraprocedural perforation",
    "bleeding": "Delayed bleeding",
    "r0": "R0 histology (0=no, 1=yes)",
    "curative_resection": "Curative resection (0=no, 1=yes)",
    "history_text": "HISTOIRE DU POLYPE: 1 MICI, 2 Radiotherapie 3 recidive post MUCO 4 divertiucles 5 appendice 6 recidive post chir",
    "history_binary": "Recidive post MUCO ou chir 1 oui 0 non",
}

PUBLIC_COLUMNS = [
    "atract",
    "study_year_index",
    "operator_id_public",
    "operator_experience",
    "age_years_topcoded",
    "sex",
    "asa",
    "anticoagulants",
    "antiplatelets",
    "location_group",
    "major_diameter_mm",
    "minor_diameter_mm",
    "surface_mm2",
    "lesion_type",
    "macronodule",
    "jnet_group",
    "conecct_group",
    "fibrosis",
    "recurrence_history",
    "procedure_duration_min",
    "speed_mm2_min",
    "r0",
    "perforation",
    "delayed_bleeding",
    "curative_resection",
]

SEX_MAP = {
    1: "male",
    2: "female",
}

EXPERIENCE_MAP = {
    1: "beginner_<50",
    2: "advanced_50_200",
    3: "expert_>200",
}

LOCATION_MAP = {
    0: "rectum",
    1: "sigmoid_left_colon",
    3: "flexures",
    4: "transverse_colon",
    5: "right_colon_cecum",
    6: "ileocecal_valve",
}

RAW_LOCATION_TO_SIMPLIFIED_MAP = {
    0: 0,
    1: 1,
    2: 1,
    3: 3,
    4: 4,
    5: 4,
    6: 3,
    7: 5,
    8: 6,
    9: 5,
}

EXCLUDED_RAW_LOCATION_CODES = {10}

LESION_TYPE_MAP = {
    1: "protruding_or_other",
    2: "lst_granular",
    3: "lst_nongranular",
    "recurrence on scar": "protruding_or_other",
}

CONECCT_MAP = {
    2: "is",
    3: "iia",
    4: "iic",
    5: "iic_plus",
}

JNET_MAP = {
    1: "jnet_i",
    2: "jnet_iia",
    3: "jnet_iib",
    4: "jnet_iii",
}

FIBROSIS_MAP = {
    0: "F0",
    1: "F1",
    2: "F2",
}

PHYSICIAN_PUBLIC_MAP = {
    "PIOCHE": "operator_01",
    "RIVORY": "operator_02",
    "LUPU": "operator_03",
    "ROSTAIN": "operator_04",
}

PUBLIC_OPERATORS = [
    "operator_01",
    "operator_02",
    "operator_03",
    "operator_04",
    "operator_other",
]

PRIMARY_OPERATOR_YEAR_MIN_PER_ARM = 3
PRIMARY_SUPPORT_SCOPE = "operator_year"
PRIMARY_PS_STRUCTURE = "operator_plus_year"
LARGE_LESION_CUTOFF_MM = 50
NN_MATCH_RATIO = 1
NN_MATCH_SCOPE = "operator_only"
PRIMARY_MATCH_CALIPER_GRID = [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60]
PRIMARY_MATCH_MAX_SMD = 0.07
PRIMARY_MATCH_MAX_OPERATOR_SHARE = 0.60

BINARY_OUTCOMES = ["r0", "perforation", "delayed_bleeding"]

CAUSAL_SPECS = {
    "primary_conecct": {
        "ps_rhs": (
            "bs(major_diameter_mm, df=4, include_intercept=False) + "
            "C(location_group) + C(lesion_type) + macronodule + C(conecct_group) + "
            "C(fibrosis) + recurrence_history"
        ),
        "balance_covariates": [
            "study_year_index",
            "operator_id_public",
            "major_diameter_mm",
            "location_group",
            "lesion_type",
            "macronodule",
            "conecct_group",
            "fibrosis",
            "recurrence_history",
        ],
        "speed_required_columns": [
            "study_year_index",
            "operator_id_public",
            "major_diameter_mm",
            "location_group",
            "lesion_type",
            "macronodule",
            "conecct_group",
            "fibrosis",
            "recurrence_history",
            "speed_mm2_min",
            "procedure_duration_min",
            "surface_mm2",
        ],
        "speed_aug_rhs": (
            "C(operator_id_public) + C(study_year_index) + "
            "bs(major_diameter_mm, df=4, include_intercept=False) + "
            "C(location_group) + C(lesion_type) + macronodule + "
            "C(conecct_group) + C(fibrosis) + recurrence_history"
        ),
        "binary_required_columns": [
            "study_year_index",
            "operator_id_public",
            "major_diameter_mm",
            "location_group",
            "lesion_type",
            "macronodule",
            "conecct_group",
            "fibrosis",
            "recurrence_history",
            "age_years_topcoded",
            "sex",
            "asa",
            "anticoagulants",
            "antiplatelets",
        ],
        "binary_aug_rhs": (
            "C(operator_id_public) + C(study_year_index) + age_years_topcoded + "
            "C(sex) + C(asa) + anticoagulants + antiplatelets + "
            "bs(major_diameter_mm, df=4, include_intercept=False) + "
            "C(location_group) + C(lesion_type) + macronodule + "
            "C(conecct_group) + C(fibrosis) + recurrence_history"
        ),
    },
}

PRIMARY_CAUSAL_SPEC = "primary_conecct"

BANNED_RAW_COLUMNS = [
    "ID",
    "NOM",
    "PRENOM",
    "Initiales Prenom",
    "Date examen",
    "Comments",
    "Raison souci technique ATRACT",
]

BANNED_STRING_TOKENS = [
    "BONNIAUD",
    "PIOCHE",
    "RIVORY",
    "ROSTAIN",
    "LUPU",
    "LAMBIN",
    "LAFEUILLE",
    "FABRITIUS",
    "LYON_",
]

DATA_DICTIONARY_ROWS = [
    ("atract", "integer", "0 = non-ATRACT traction, 1 = ATRACT traction", "ATRACT (0= non, 1= oui)"),
    ("study_year_index", "integer", "Ordered study year index; 1 is earliest year in the cohort", "Date examen"),
    ("operator_id_public", "string", "Pseudonymized operator identifier; rare operators collapsed to operator_other", "Physician"),
    ("operator_experience", "string", "Operator experience group from the source workbook", "Operateur simplifie"),
    ("age_years_topcoded", "integer", "Age in years, top-coded at 90", "Âge"),
    ("sex", "string", "Biological sex category from the source workbook", "Gender (1 male, 2 woman)"),
    ("asa", "integer", "ASA physical status class", "ASA"),
    ("anticoagulants", "integer", "0 = no, 1 = yes", "Anticoaggulants"),
    ("antiplatelets", "integer", "0 = no, 1 = yes", "Antiaggrégants"),
    ("location_group", "string", "Six-level lesion location grouping used in the public analysis", "Localisation simplifiee"),
    ("major_diameter_mm", "float", "Largest specimen diameter in millimeters", "Major diameter (mm)"),
    ("minor_diameter_mm", "float", "Minor specimen diameter in millimeters", "Minor Diameter (mm)"),
    ("surface_mm2", "float", "Dissected specimen surface area in square millimeters", "Surface"),
    ("lesion_type", "string", "Normalized lesion morphology category", "Type de lésion (1 Polype, 2 LST Granulaire, 3 LST non granulaire)"),
    ("macronodule", "integer", "0 = absent, 1 = present", "Macronodule"),
    ("jnet_group", "string", "Normalized JNET class used in sensitivity analyses", "JNET (1= JNET  I, 2= JNET IIa, 3= JNET IIB, 4= JNET III)"),
    ("conecct_group", "string", "Normalized CONECCT class used in modeling", "CONECCT (1=1H, 2=Is, 3=IIA, 4=IIC, 5=IIC+, 6=III, 7=OE)"),
    ("fibrosis", "string", "Normalized fibrosis category", "Fibrosis"),
    ("recurrence_history", "integer", "0 = no recurrence/scar history, 1 = recurrence or scar history", "HISTOIRE DU POLYPE... / Recidive post MUCO ou chir 1 oui 0 non"),
    ("procedure_duration_min", "float", "Procedure duration in minutes", "Duration of procedure Procédure (min)"),
    ("speed_mm2_min", "float", "Dissection speed in square millimeters per minute", "Vitesse"),
    ("r0", "integer", "0 = non-R0, 1 = R0 histology", "R0 histology (0=no, 1=yes)"),
    ("perforation", "integer", "0 = no intraprocedural perforation, 1 = perforation", "Intraprocedural perforation"),
    ("delayed_bleeding", "integer", "0 = no delayed bleeding, 1 = delayed bleeding", "Delayed bleeding"),
    ("curative_resection", "integer", "0 = non-curative resection, 1 = curative resection", "Curative resection (0=no, 1=yes)"),
]
