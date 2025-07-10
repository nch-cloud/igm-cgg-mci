import glob, json, os, copy
import argparse

## Hard-coded Reference Paths

methy_ref_file_v11 = os.path.join(os.path.dirname(os.path.realpath(__file__)),"methyl_v11.csv")
methyl_ref = {}

def prep_methyl_ref_v11():
    import pandas as pd
    v11_ref = pd.read_csv(methy_ref_file_v11)
    class_key = v11_ref["class_code"]
    class_value = v11_ref["class_string"]
    family_key = v11_ref["family_code"]
    family_value = v11_ref["family_string"]
    v11_dict = {}
    for i in range(len(class_key)):
        adj_key = f"{class_key[i]}".replace("/"," ").replace(","," ").replace("  "," ").strip().replace(" ","_").replace("MCF_","").replace("MTF_","").replace("MTGF_","")
        while "__" in adj_key:
            adj_key = adj_key.replace("__","_")
        v11_dict[adj_key]=f"{class_value[i]}".replace("methylation class ","").strip()
    for i in range(len(family_key)):
        adj_key = f"{family_key[i]}".replace("/"," ").replace(","," ").replace("  "," ").strip().replace(" ","_").replace("MCF_","").replace("MTF_","").replace("MTGF_","")
        while "__" in adj_key:
            adj_key = adj_key.replace("__","_")
        if len(adj_key) > 0:
            v11_dict[adj_key]=f"{family_value[i]}".replace("methylation class family ","").strip()
    methyl_ref["v11"] = v11_dict

## JSON Deduplication Functions

def handle_duplicates(pairs):
    result = {}
    for key, value in pairs:
        show = False
        iter = 0
        new_key = key
        while new_key in result:
            iter += 1
            new_key = f"{key}_{iter}"
            show = True
        result[new_key] = value
        if show:
            print(f"Deduplicated {key} to {new_key}")
    return result

## Prep Methods

def get_dir_jsons(target_dir:str, log:bool=False):
    #print(target_dir)
    if log:
        print(f"Globbing JSONs in {target_dir}...")
    # Get all JSON files in the target directory and return them as a list.
    jsons = glob.glob(os.path.join(target_dir,"*.json"))

    if log:
        print(f"Found {len(jsons)} files.")

    return jsons

def sort_jsons(json_dir_list:list, log:bool = False):
    # Restructures JSONs into a dictionary keyed by samples.
    # Values are copies of the 'blank dict' below, which is populated with JSON data.
    blank_dict = {"cog":None, "tumor_normal":None,  "methyl_igm":None, "methyl_v11":None, "methyl_v12":None, "archer_fusion":None, "methyl_v11_raw":None, "methyl_v12_raw":None, "methyl_igm_raw":None}
    json_dicts = {}

    json_type = None
    
    json_sizes = {}

    json_list = []
    for dir in json_dir_list:
        #print(dir)
        json_list.extend(get_dir_jsons(dir))

    for i in json_list:
        json_type = None

        try:
            with open(i, 'r') as json_file:
                json_data = json.load(json_file, object_pairs_hook=handle_duplicates)
            if 'subject_id' in json_data:
                subject = json_data['subject_id']
                report_type = json_data['report_type']
                if 'meth' in report_type.lower():
                    if 'report_version' in json_data:
                        if 'IGM' in json_data['report_version']:
                            json_type = "methyl_igm"
                        elif 'v12' in json_data['report_version']:
                            json_type = "methyl_v12"
                        elif 'v11' in json_data['report_version']:
                            json_type = "methyl_v11"
                elif report_type == "archer_fusion":
                    json_type = "archer_fusion"
                elif report_type == "tumor_normal":
                    json_type = "tumor_normal"
                else:
                    if log:
                        print(f"Unknown report type: {report_type} [{i}]. Skipping...")
            elif 'upi' in json_data:
                # COG File
                subject = json_data['upi']
                json_type = "cog"
            elif 'meta_data' in json_data and 'report_title' in json_data['meta_data'] and "Methylation" in json_data['meta_data']['report_title']:
                # Methylation data
                subject = os.path.basename(i).replace("-","_").split("_")[1]
                if "v12" in json_data['meta_data']['report_title']:
                    json_type = "methyl_v12_raw"
                elif "IGM" in json_data['meta_data']['report_title']:
                    json_type = "methyl_igm_raw"
                else:
                    json_type = "methyl_v11_raw"
            else:
                if log:
                    print(f"Skipping {i}...")
                continue

            if json_type is not None and json_type in blank_dict:
                if log:
                    print(f"Processing {i} as {json_type}.")
                if subject not in json_dicts:
                    json_dicts[subject] = copy.copy(blank_dict)
                    json_sizes[subject]={}
                json_size = os.path.getsize(i)
                if json_dicts[subject][json_type] is not None:
                    if json_size == json_sizes[subject][json_type]:
                        if log:
                            print(f"{subject} {json_type} is duplicated. Skipping overwrite.")
                        continue
                    elif json_size < json_sizes[subject][json_type]:
                        if log:
                            print(f"New {subject} {json_type} is smaller. Skipping overwrite.")
                        continue
                    else:
                        if log:
                            print(f"Overwriting {subject} {json_type} with newer version.")
            

                json_sizes[subject][json_type]=json_size
                json_dicts[subject][json_type]=json_data
        except OSError as e:
            print(f"ERROR: {str(e)}. File {i}'s data will be missing from outputs.")

    return json_dicts

# Clean-Up Methods

def replace_blank_fields(data:dict, blank_field_placeholder:str = "."):
    # Goes through and replaces values that are an empty string with a pre-specified value.
    # The purpose of this is to show which fields are deliberately left blank, vs. which are due to missing data files. 
    for i in data:
        if type(data[i]) is str and data[i] == "":
            data[i] = blank_field_placeholder
    return data

def standardize_variant_notation(data:dict,debug:bool=False,
                                fields:list=["TN_Germline_Path","TN_Germline_LikelyPath","TN_Germline_VUS","TN_Somatic_Tier1","TN_Somatic_Tier2","TN_Somatic_Tier3"]
                                ):
    # Goes through and gets the longest-form version of each variant notation, and then expands all matching variants to conform to it
    if debug:
        print("\nSynonymizing variant notation for fields:\n\t")
        print("\n\t".join(fields))

    # Get all variants
    variants = set()
    for i in data:
        for j in fields:
            if j in i:
                if len(i[j]) <= 1:
                    pass
                else:
                    for v in i[j].split(";"):
                        variants.add(v)
    
    # Split into first three fields (A) and last field (B)
    parts_A = []
    parts_B = []

    for i in variants:
        parts = i.split(" ")
        last_part = parts.pop(-1)
        parts_A.append(" ".join(parts[0:3]))
        parts_B.append(last_part)

    convert_dict_A = {}

    # If a field is a longer version of another field, make one into a conversion of the other
    for i in range(len(parts_A)):
        for j in range(len(parts_A)):
            if parts_A[j].startswith(parts_A[i]) and len(parts_A[j]) == len(parts_A[i]):
                convert_dict_A[parts_A[i]] = parts_A[j]
                parts_A[i] = parts_A[j]

    expand_dict = {}

    # Make a dictionary to replace the original part B.
    for i in range(len(parts_A)):
        if parts_A[i] not in expand_dict:
            expand_dict[parts_A[i]]={}
        expand_dict[parts_A[i]][len(parts_B[i])]=parts_B[i]
    
    convert_dict_B={}

    # Populate dict with the longest part B for each option
    for i in range(len(parts_A)):
        max_len = max(expand_dict[parts_A[i]].keys())
        parts_B[i] = expand_dict[parts_A[i]][max_len]
        convert_dict_B[parts_A[i]]=parts_B[i]

    # Apply updated data
    for i in range(len(data)):
        for j in fields:
            if j in data[i]:
                if len(data[i][j]) <= 1:
                    pass
                else:
                    variants = []
                    for v in data[i][j].split(";"):
                        var_str = " ".join(v.split(" ")[0:3])
                        update_A = convert_dict_A[var_str]
                        if (update_A) in convert_dict_A:
                            update_A = convert_dict_A[update_A]
                        update_B = f"{update_A} {convert_dict_B[update_A]}"
                        variants.append(update_B)
                        if (v != update_B and debug):
                            print(f"{v} -> {update_B}")
                    data[i][j]=";".join(variants)

    return data

def methyl_to_ref(methyl_raw:str):
    for x in [".",",",";",":","-","(",")","[","]","/","&","and"]:
        methyl_ref = methyl_raw.replace(x, "_")
    while "__" in methyl_ref:
        methyl_ref = methyl_ref.replace("__","_")
    methyl_ref = "".join(sorted(methyl_ref.lower().split("_")))
    return methyl_ref

def standardize_methylation_class(data:dict,debug:bool=False,
    fields:list=["Methylation_Superfamily","Methylation_Family","Methylation_Class","Methylation_Subclass"]
    ):
    # Goes through and determines a standard version of methylation class
    if debug:
        print("\nSynonymizing methylation notation for fields:")
        print('\t' + "\n\t".join(fields))

    methyl_classes = {}

    for i in range(len(data)):
        if data[i] is None:
            continue
        for j in fields:
            if j is None or j == "":
                continue
            if j in data[i]:
                methyl_data = data[i][j][::-1].strip()[::-1].strip().replace("haem","hem").replace("Haem","Hem").replace("paed","ped").replace("Paed","ped").replace("_"," ")
                while len(methyl_data) > 0 and methyl_data[-1] in [".",",",";",":"," "]:
                    methyl_data = methyl_data[0:-1]

                #print(f'{methyl_data} -> {data[i][j]}')
                
                data[i][j] = methyl_data

                methyl_ref = methyl_to_ref(methyl_data)
                if methyl_ref not in methyl_classes:
                    methyl_classes[methyl_ref]={}
                if methyl_data not in methyl_classes[methyl_ref]:
                    methyl_classes[methyl_ref][methyl_data]=1
                else:
                    methyl_classes[methyl_ref][methyl_data]+=1

    methyl_convert = {}

    for i in methyl_classes:
        print(i)
        methyl_max = max(list(methyl_classes[i].values()))
        for j in methyl_classes[i]:
            print(f"\t{methyl_classes[i][j]} : {j}")
            if methyl_classes[i][j] == methyl_max:
                methyl_convert[i]=j

    for i in range(len(data)):
        if data[i] is None:
            continue
        for j in fields:
            if j in data[i]:
                new_methyl = methyl_convert[methyl_to_ref(data[i][j])]
                if len(new_methyl) > 1:
                    new_methyl = new_methyl[0].upper() + new_methyl[1:]
                if (data[i][j] != new_methyl):
                    if debug :
                        print(f"{data[i][j]} -> {new_methyl}")
                    data[i][j] = new_methyl

    return data

# Clinical JSONs processing

##COG JSON

def parse_cog_json(cog_json, out_dict:dict={}):
    # Takes in the sample JSONs dictionary and processes the COG JSON if it is present.
    if cog_json is None:
        return out_dict
    forms_dict = {}
    follow_ups = []
    for i in cog_json['forms']:
        form_name = i['form_id']
        if form_name == "FOLLOW_UP":
            # Follow-up can have multiple forms, which all have key 'data' in the COG JSON.
            fup_keys = i.keys()
            for f in fup_keys:
                if f.startswith("data"):
                    fup_dict = {}
                    form_data = i[f]
                    for j in form_data:
                        field_id = j['form_field_id']
                        fup_dict[field_id]=j
                    follow_ups.append(fup_dict)
        else:
            forms_dict[form_name]={}
            form_data = i['data']
            for j in form_data:
                field_id = j['form_field_id']
                forms_dict[form_name][field_id]=j
    # If there were any follow-up forms, replace the single form in the forms-dict with the collection of follow-ups.
    if len(follow_ups) > 0:
        forms_dict["FOLLOW_UP"]=follow_ups

    # Demography form
    if "DEMOGRAPHY" in forms_dict:
        out_dict["Birth_Date"]=forms_dict["DEMOGRAPHY"]["DM_BRTHDAT"]["value"]
        out_dict["Ethnicity"]=forms_dict["DEMOGRAPHY"]["DM_ETHNIC"]["value"]
        out_dict["Race"]=forms_dict["DEMOGRAPHY"]["DM_CRACE"]["value"]
        out_dict["Sex"]=forms_dict["DEMOGRAPHY"]["DM_SEX"]["value"]
        out_dict["Country_of_Residence"]=forms_dict["DEMOGRAPHY"]["SC_SCORRES_CNTRYRES"]["value"]
    # Diagnosis form
    if "COG_UPR_DX" in forms_dict:
        out_dict["Diagnosis_ID"]=forms_dict["COG_UPR_DX"]["ADM_DX_CD_SEQ"]["value"]
        out_dict["Enrolled_Dx"]=forms_dict["COG_UPR_DX"]["PRM_TU_DX_TXT"]["value"]
        out_dict["Date_of_Diagnosis"]=forms_dict["COG_UPR_DX"]["DX_DT"]["value"]
        out_dict["Primary_Site_Code"]=forms_dict["COG_UPR_DX"]["TOPO_ICDO"]["value"]
        out_dict["Primary_Site_Term"]=forms_dict["COG_UPR_DX"]["TOPO_TEXT"]["value"]
        out_dict["Initial_Dx_Code"]=forms_dict["COG_UPR_DX"]["MORPHO_ICDO"]["value"]
        out_dict["Initial_Dx_Term"]=forms_dict["COG_UPR_DX"]["MORPHO_TEXT"]["value"]
        out_dict["Registry_Stage_Code"]=forms_dict["COG_UPR_DX"]["REG_STAGE_CODE_TEXT"]["value"]
    # Registry data form
    if "REGISTRY_DATA" in forms_dict:
        out_dict["Date_of_Death"]=forms_dict["REGISTRY_DATA"]["DEATH_DOC_DATE"]["value"]
    # Final Diagnosis form
    if "FINAL_DIAGNOSIS" in forms_dict:
        out_dict["Dx_Morpho_Code"]=forms_dict["FINAL_DIAGNOSIS"]["PRM_CA_DX_ICD_O_CD"]["value"]
        out_dict["Primary_Dx_Disease_Group"]=forms_dict["FINAL_DIAGNOSIS"]["PRIMDXDSCAT"]["value"]
    # Treatment confirmation form
    if "TREATMENT_CONFIRMATION" in forms_dict:
        out_dict["Enrolled_on_Prev_COG_Study"]=forms_dict["TREATMENT_CONFIRMATION"]["PT_OTH_ENROLLM_IND_2"]["value"]
    # On-Study Diagnosis (CNS) form
    if "ON_STUDY_DX_CNS" in forms_dict:
        out_dict["Tumor_Grade"]=forms_dict["ON_STUDY_DX_CNS"]["TUMOR_GP_ST"]["value"]
        out_dict["Tumor_M_Stage"]=forms_dict["ON_STUDY_DX_CNS"]["CNSTMRMSTG"]["value"]
        out_dict["Cerebrospinal_Fluid_Status"]=forms_dict["ON_STUDY_DX_CNS"]["CSFCYTLGY"]["value"]
        out_dict["Spine_at_Diagnosis"]=forms_dict["ON_STUDY_DX_CNS"]["CNSSPNDXSTATUS"]["value"]
        out_dict["Had_Surgical_Resection"]=";".join([forms_dict["ON_STUDY_DX_CNS"]["SURGBXRCTPFM"]["value"],
            forms_dict["ON_STUDY_DX_CNS"]["TUM_RES_EXT_TP"]["value"],forms_dict["ON_STUDY_DX_CNS"]["OTX_SURG_RESECT_TXT"]["value"]])
        out_dict["Residual_Tumor"]=forms_dict["ON_STUDY_DX_CNS"]["RESI_MALI_POST_SURG_MEAS"]["value"]
    # NCI_MCI_FUP (?) form
    if "NCI_MCI_FUP" in forms_dict:
        out_dict["Has_Molecular_Reports"]=forms_dict["NCI_MCI_FUP"]["MCRPTRCVD"]["value"]
        out_dict["Trial_Enrolled_Using_Results"]=forms_dict["NCI_MCI_FUP"]["PTNTENRLSEQELIGTREATASGNIND"]["value"]
        out_dict["Therapy_Matched_By_Sequencing"]=forms_dict["NCI_MCI_FUP"]["PTNTMOLSEQVARINDMCHTXTRLENR"]["value"]
        out_dict["Dx_Refined_by_Testing"]=forms_dict["NCI_MCI_FUP"]["FNLDXMOLANLSUPDOTCM"]["value"]
    # Follow-Up form
    if "FOLLOW_UP" in forms_dict:
        for fup_form in forms_dict["FOLLOW_UP"]:
            append_string(fup_form["REP_EVAL_PD_TP"]["value"], "APEC14B1_Reporting_Period", out_dict)
            if fup_form["PT_INF_CU_FU_COL_IND"]["value"].lower() == "yes":
                fup_obt = f"{fup_form["PT_INF_CU_FU_COL_IND"]["value"]} ({fup_form["PT_FU_BEGDT"]["value"]}-{fup_form["PT_FU_END_DT"]["value"]})"
            else:
                fup_obt = fup_form["PT_INF_CU_FU_COL_IND"]["value"]
            append_string(fup_obt, "FollowUp_Obtained_for_Period", out_dict)
            append_string(fup_form["PT_VST"]["value"], "Vital_status", out_dict)
            frontline_treatments=";".join(get_frontline_treatments([fup_form["FSTLNTXINIDXADM"]["value"],
                fup_form["FSTLNTXINIDXADMCAT_A1"]["value"],
                fup_form["FSTLNTXINIDXADMCAT_A2"]["value"],
                fup_form["FSTLNTXINIDXADMCAT_A3"]["value"],
                fup_form["FSTLNTXINIDXADMCAT_A4"]["value"],
                fup_form["FSTLNTXINIDXADMCAT_A5"]["value"],
                fup_form["FSTLNTXINIDXADMCAT_A6"]["value"],
                fup_form["FSTLNTXINIDXADMOS"]["value"]]))
            if len(frontline_treatments) > 0 or "Frontline_Treatment_Received" not in out_dict:
                append_string(frontline_treatments, "Frontline_Treatment_Received", out_dict)
            append_string(fup_form["DZ_EXM_REP_IND_2"]["value"], "Disease_Status_Evaluated_During_Interval", out_dict)
            append_string(fup_form["COMP_RESP_CONF_IND_3"]["value"],"Achieved_Complete_Remission", out_dict)
            append_string(fup_form["DZ_REL_PROG_IND3"]["value"],"Developed_First_Relapse_or_Progression", out_dict)
            append_string(fup_form["NEW_CA_DX_IND_3"]["value"], "Dx_New_Primary_or_MDS", out_dict)
            append_string(fup_form["PT_FU_ANNIV_REACH_IND"]["value"], "Patient_Reached_Tenth_Anniv", out_dict)
            append_string(fup_form["PT_LOST_FU_IND_2"]["value"], "Confirmed_Lost_to_FollowUp", out_dict)
            append_string(fup_form["PT_FOL_CON_IND"]["value"], "Plans_To_Continue_Tracking_Outcome", out_dict)
            append_string(fup_form["PTWDRWCSNTFUENDRPDIND"]["value"], "Withdrew_APEC14B1_Consent", out_dict)
    # On-study diagnosis (STS) form
    if "ON_STUDY_DX_SOFT_TISSUE_SARCOMA" in forms_dict:
        out_dict["Procedure_Type"]=forms_dict["ON_STUDY_DX_SOFT_TISSUE_SARCOMA"]["SURG_RESECT_EXT_TP"]["value"]
        if out_dict["Procedure_Type"] == "Other":
            out_dict["Procedure_Type"] = f'{out_dict["Procedure_Type"]};{forms_dict["ON_STUDY_DX_SOFT_TISSUE_SARCOMA"]["SURG_PROC_O_SPEC_TXT"]["value"]}'
    # CNS Chemo Treatment form
    if "TX_CHEMO_CNS" in forms_dict:
        out_dict["Treated_but_not_Enrolled"]=forms_dict["TX_CHEMO_CNS"]["TX_RCVD_YES_NO"]["value"]
        out_dict["COG_Anti_Cancer_Treatment"]=forms_dict["TX_CHEMO_CNS"]["COG_ID_ENUM"]["value"]
        if out_dict["COG_Anti_Cancer_Treatment"].lower() == "other":
            out_dict["COG_Anti_Cancer_Treatment"] = ";".join([out_dict["COG_Anti_Cancer_Treatment"],forms_dict["TX_CHEMO_CNS"]["COG_ID_OTHER"]["value"],forms_dict["TX_CHEMO_CNS"]["PRI_TX_RGM_SPEC"]["value"]])
        out_dict["Non_COG_Anti_Cancer_Treatment"]=forms_dict["TX_CHEMO_CNS"]["NPROT_TX_ADM_IND_3"]["value"]
        if out_dict["Non_COG_Anti_Cancer_Treatment"].lower() == "other":
            out_dict["Non_COG_Anti_Cancer_Treatment"] = ";".join([out_dict["Non_COG_Anti_Cancer_Treatment"],forms_dict["TX_CHEMO_CNS"]["NPROT_TX_ADM_NM"]["value"],forms_dict["TX_CHEMO_CNS"]["NPROT_TX_ADM_SPEC"]["value"]])
        out_dict["Chemotherapy"]=";".join(get_chemo_drugs(forms_dict["TX_CHEMO_CNS"])) # Uses a helper function to join all these
    # Death form
    if "DEATH" in forms_dict:
        out_dict["Primary_Cause_of_Death"]=forms_dict["DEATH"]["PT_DEATH_PRM_RSN"]["value"]
    # Radiation therapy form
    if "RADIATION_THERAPY" in forms_dict:
        out_dict["Radiation_Therapy"]=";".join(get_radiation_types(forms_dict["RADIATION_THERAPY"])) # Uses a helper function to join all these
    # Relapse/Prognosis (CNS) form
    if "RLP_PROG_CNS" in forms_dict:
        out_dict["Relapse_Status"]=forms_dict["RLP_PROG_CNS"]["PROG_REL_STAT"]["value"]
        out_dict["Relapse_Date"]=forms_dict["RLP_PROG_CNS"]["DZ_RECUR_PROG_DX_DT"]["value"]
        out_dict["Relapse_Site"]=forms_dict["RLP_PROG_CNS"]["MET_REL_PROG_LOC_CATE_A1"]["value"]
    if "CNS_DIAGNOSIS_DETAIL" in forms_dict:
        out_dict["CNS_Diagnosis_Category"]=forms_dict["CNS_DIAGNOSIS_DETAIL"]["MH_MHCAT_CNSDXCAT"]["value"]
        integrated_dx = [] # This is a list to handle possible edge-cases where more than one integrated diagnosis field has data
        for field in forms_dict["CNS_DIAGNOSIS_DETAIL"]:
            if field.startswith("MH_MHSCAT_CNSDXINTGRT"):
                if len(forms_dict["CNS_DIAGNOSIS_DETAIL"][field]["value"]) > 1:
                    integrated_dx.append(forms_dict["CNS_DIAGNOSIS_DETAIL"][field]["value"])
        out_dict["CNS_Integrated_Diagnosis"]=";".join(integrated_dx)

    return out_dict

def get_chemo_drugs(cns_chemo_dict:dict):
    # Helper function to concatenate chemo drug fields
    drugs = []
    for i in cns_chemo_dict:
        if i.startswith('AGT_ADM_NM'):
            if cns_chemo_dict[i]['value'] == 'checked':
                drugs.append(cns_chemo_dict[i]['SASLabel'])
    return drugs

def get_radiation_types(cns_rad_dict:dict):
    # Helper function to concatenate radiation treatment fields
    rad = []
    for i in cns_rad_dict:
        if i.startswith('RT_TX_TP'):
            if cns_rad_dict[i]['value'] == 'checked':
                rad.append(cns_rad_dict[i]['SASLabel'])
    return rad

def get_frontline_treatments(frontline_treatments:list):
    # Helper function to convert frontline treatments to list of strings
    treatments_ref = ["Chemotherapy/Immunotherapy", "Radiation Therapy", "Stem Cell Transplant", "Surgery", "Cellular Therapy"]
    treatments_out = []
    for i in range(len(treatments_ref)):
        if frontline_treatments[i+1].lower() == "checked":
            treatments_out.append(treatments_ref[i])
    if frontline_treatments[7].lower() == "checked":
        treatments_out.append(frontline_treatments[8])
    if len(treatments_out) == 0:
        treatments_out = [frontline_treatments[0]]
    return treatments_out

def append_string(new_data:str, data_key:str, data_dict:dict, separator:str=";", blank_field_indicator:str=""):
    if new_data is None or new_data == "":
        new_data = blank_field_indicator
    if data_key in data_dict:
        data_dict[data_key] = f"{data_dict[data_key]}{separator}{new_data}"
    else:
        data_dict[data_key] = new_data

## IGM Molecular

def parse_molecular_generic(generic_json, out_dict:dict={}):
    if generic_json is None:
        return out_dict
    if "Cellularity" not in out_dict:
        out_dict["Cellularity"] = generic_json["percent_tumor"]
        out_dict["Necrosis"] = generic_json["percent_necrosis"]
    if "Disease_Group" not in out_dict:
        out_dict["Disease_Group"]=generic_json["disease_group"]
    if "Indication" not in out_dict:
        out_dict["Indication"]=generic_json["indication_for_study"]
    return out_dict

### Tumor-Normal Exome

def parse_tumor_normal_json(tn_json, out_dict:dict={}):
    if tn_json is None:
        return out_dict
    
    if 'TN_Version' in out_dict and out_dict['TN_Version'] is not None:
        print(f"WARNING! Duplicate TN JSONs present for {out_dict['Sample']}")

    out_dict['TN_Version'] = tn_json['version']
    out_dict['TN_Germline_Result']='Negative'
    out_dict['TN_Somatic_Result']='Negative'
    out_dict['TN_Germline_CNV_Result']='Negative'
    out_dict['TN_Somatic_CNV_Result']='Negative'
                
    variants={'TN_Germline_Path':[],'TN_Germline_LikelyPath':[],'TN_Germline_VUS':[],
              'TN_Germline_CNV_Tier1-2':[],'TN_Germline_CNV_Tier3':[],
              'TN_Somatic_Tier1':[],'TN_Somatic_Tier2':[],'TN_Somatic_Tier3':[],
              'TN_Somatic_CNV_Tier1-2':[],'TN_Somatic_CNV_Tier3':[]}
    gene_changes={'TN_Germline_CNV_Gene_Loss':[],'TN_Germline_CNV_Gene_BiallelicLoss':[],
                  'TN_Germline_CNV_Gene_Gain':[],'TN_Germline_CNV_Gene_Amplification':[],
                  'TN_Somatic_CNV_Gene_Loss':[],'TN_Somatic_CNV_Gene_BiallelicLoss':[],
                  'TN_Somatic_CNV_Gene_Gain':[],'TN_Somatic_CNV_Gene_Amplification':[],
                  'TN_Germline_CNV_Gene_LOH':[], 'TN_Somatic_CNV_Gene_LOH':[],
                  'TN_Germline_CNV_Blurb':"", 'TN_Somatic_CNV_Blurb':""
                  }
    if 'somatic_results' in tn_json:
        if 'variants' in tn_json['somatic_results']:
            for v in tn_json['somatic_results']['variants']:
                var_str, tier = var_to_string(v)
                if tier == '1':
                    variants['TN_Somatic_Tier1'].append(var_str)
                elif tier == '2':
                    variants['TN_Somatic_Tier2'].append(var_str)
                else:
                    variants['TN_Somatic_Tier3'].append(var_str)
                if len(var_str) > 0:
                    out_dict['TN_Somatic_Result']='Positive'
    if 'germline_results' in tn_json:
        if 'variants' in tn_json['germline_results']:
            for v in tn_json['germline_results']['variants']:
                var_str, tier = var_to_string(v)
                if tier == 'LikelyPath':
                    variants['TN_Germline_LikelyPath'].append(var_str)
                elif tier == 'Path':
                    variants['TN_Germline_Path'].append(var_str)
                else:
                    variants['TN_Germline_VUS'].append(var_str)
                if len(var_str) > 0:
                    out_dict['TN_Germline_Result']='Positive'

    if 'somatic_cnv_results' in tn_json:
        if 'variants' in tn_json['somatic_cnv_results']:
            for v in tn_json['somatic_cnv_results']['variants']:
                var_str, tier = cnv_to_string(v)
                if tier == '1' or 'tier' == '2':
                    variants['TN_Somatic_CNV_Tier1-2'].append(var_str)
                else:
                    variants['TN_Somatic_CNV_Tier3'].append(var_str)
                out_dict['TN_Somatic_CNV_Result']='Positive'
                copy_type = v['copy_number_type'].lower()
                if 'bialle' in copy_type or 'complet' in copy_type or 'total' in copy_type:
                    gene_changes['TN_Somatic_CNV_Gene_BiallelicLoss'].extend(v['disease_associated_gene_content'])
                elif 'loh' in copy_type or 'hetero' in copy_type or 'roh' in copy_type or 'homo' in copy_type:
                    gene_changes['TN_Somatic_CNV_Gene_LOH'].extend(v['disease_associated_gene_content'])   
                elif 'gain' in copy_type:
                    gene_changes['TN_Somatic_CNV_Gene_Gain'].extend(v['disease_associated_gene_content'])
                elif 'ampli' in copy_type:
                    gene_changes['TN_Somatic_CNV_Gene_Amplification'].extend(v['disease_associated_gene_content'])
                elif 'loss' in copy_type:
                    gene_changes['TN_Somatic_CNV_Gene_Loss'].extend(v['disease_associated_gene_content'])
                else:
                    #print(f'{tn_json["subject_id"]} Invalid somatic CNV-LOH event type: {copy_type}')
                    pass
        if 'summary' in tn_json['somatic_cnv_results']:
            gene_changes["TN_Somatic_CNV_Blurb"] = ("\n".join(tn_json['somatic_cnv_results']['summary']))
            if out_dict['TN_Somatic_CNV_Result'] == 'Negative' and gene_changes["TN_Somatic_CNV_Blurb"] is not None:
                bt = gene_changes["TN_Somatic_CNV_Blurb"].strip().lower()
                if len(bt) > 0 and not any(["none detected" in bt, "not sufficient" in bt, "insufficient" in bt, "did not meet" in bt]):
                    out_dict['TN_Somatic_CNV_Result']='Positive'

    if 'germline_cnv_results' in tn_json:
        if 'variants' in tn_json['germline_cnv_results']:
            for v in tn_json['germline_cnv_results']['variants']:
                var_str, tier = cnv_to_string(v)
                if tier == '1' or 'tier' == '2':
                    variants['TN_Germline_CNV_Tier1-2'].append(var_str)
                else:
                    variants['TN_Germline_CNV_Tier3'].append(var_str)
                out_dict['TN_Germline_CNV_Result']='Positive'
                copy_type = v['copy_number_type'].lower()
                if 'bialle' in copy_type or 'complet' in copy_type or 'total' in copy_type:
                    gene_changes['TN_Germline_CNV_Gene_BiallelicLoss'].extend(v['disease_associated_gene_content'])
                elif 'loh' in copy_type or 'hetero' in copy_type or 'roh' in copy_type or 'homo' in copy_type:
                    gene_changes['TN_Germline_CNV_Gene_LOH'].extend(v['disease_associated_gene_content'])   
                elif 'gain' in copy_type:
                    gene_changes['TN_Germline_CNV_Gene_Gain'].extend(v['disease_associated_gene_content'])
                elif 'ampli' in copy_type:
                    gene_changes['TN_Germline_CNV_Gene_Amplification'].extend(v['disease_associated_gene_content'])
                elif 'loss' in copy_type or 'del' in copy_type:
                    gene_changes['TN_Germline_CNV_Gene_Loss'].extend(v['disease_associated_gene_content'])
                else:
                    #print(f'{tn_json["subject_id"]} Invalid germline CNV-LOH event type: {copy_type}')
                    pass
        if 'summary' in tn_json['germline_cnv_results']:
            gene_changes["TN_Germline_CNV_Blurb"] = ("\n".join(tn_json['germline_cnv_results']['summary']))
            if out_dict['TN_Germline_CNV_Result'] == 'Negative' and gene_changes["TN_Germline_CNV_Blurb"] is not None:
                bt = gene_changes["TN_Germline_CNV_Blurb"].strip().lower()
                if len(bt) > 0 and not any(["none detected" in bt, "not sufficient" in bt, "insufficient" in bt, "did not meet" in bt]):
                    out_dict['TN_Germline_CNV_Result']='Positive'
    
    for i in variants:
        #for j in range(len(variants[i])):
        #    variant = variants[i][j]
        #    if "(" in variant or ")" in variant:
        #        print(variant)
        out_dict[i]=";".join(variants[i])

    for i in gene_changes:
        if "Blurb" in i:
            gene_string = gene_changes[i]
        else:
            gene_set = set(gene_changes[i])
            gene_set.discard("N/A")
            gene_set.discard("NA")
            gene_set.discard("n/a")
            gene_set.discard("")
            gene_string = ";".join(sorted(list(gene_set)))
        out_dict[i] = gene_string

    out_dict = parse_molecular_generic(tn_json, out_dict)
    return out_dict

def get_tier(tier_str:str):
    # Process various strings to tiering
    tier_str = tier_str.lower()
    if 'tier 3' in tier_str or 'tier iii' in tier_str:
        tier = '3'
    elif 'tier 2' in tier_str or 'tier ii' in tier_str:
        tier = '2'
    elif 'tier 1' in tier_str or 'tier i' in tier_str:
        tier = '1'
    elif 'likely' in tier_str:
        tier = 'LikelyPath'
    elif 'path' in tier_str:
        tier = 'Path'
    else:
        tier = 'VUS'
    return tier

def var_to_string(info:dict):
    # Process variant notation into string
    gene = info['gene'].replace(" ","")
    trans = info['transcript'].replace(" ","")
    nuc_change = info['nucleotide_change'].replace(" ","")
    if 'predicted_protein_change' in info:
        prot_change = info['predicted_protein_change'].replace(" ","")
    else:
        prot_change = 'p.?'
    if '=' in prot_change:
        amino = prot_change.split(".")[1][0:3]
        prot_change = prot_change.replace("=", amino)
    if "*" in prot_change:
        prot_change = prot_change.replace("*","Ter")
    var_str = f'{gene} {trans} {nuc_change} {prot_change}'.replace("(","").replace(")","").replace(";",",")
    tier = get_tier(info['interpretation']['value'])
    return var_str, tier

def cnv_to_string(info:dict):
    cnv_type = info['copy_number_type']
    if info['genomic_change']['chromosome'] is None: #This is most likely a special type of change
        var_str = info['cytogenetic_locus'].title().replace("Near ","Near-")
        #print(var_str)
    else:
        if cnv_type.lower().startswith('whole chrom'):
            var_str = f"{info['genomic_change']['chromosome']} {cnv_type}"
        else:
            cnv_region = f"{info['genomic_change']['chromosome']}:{info['genomic_change']['start']}-{info['genomic_change']['end']}"
            var_str = f"{cnv_region} {cnv_type}"
            if 'focal' in cnv_type.lower():
                genes = [info['disease_associated_gene_content'][0]]
                if '(exon' in cnv_type.lower():
                    var_str = var_str.replace('(exon',f"({genes[0]} exon")
                else:
                    if len(info['disease_associated_gene_content']) > 1:
                        genes.append(info['disease_associated_gene_content'][-1])
                    var_str = f"{var_str} ({','.join(genes)})"
        var_str = var_str.replace("\\n","").replace("\t"," ").replace("Loss of Heterozygosity","LOH").replace("Loss-of_Heterozygosity","LOH").replace("(LOH)","LOH").replace("LOH LOH","LOH")
        while "  " in var_str:
            var_str = var_str.replace("  ", " ")
    tier = get_tier(info['interpretation']['value'])
    return var_str, tier

### Methylation

def parse_methyl_json(methyl_json, type:str=None, out_dict:dict={}):
    # Parses IGM methylation outputs
    if methyl_json is None:
        return out_dict
    out_dict = parse_molecular_generic(methyl_json, out_dict)
    out_dict['Methylation_Version'] = methyl_json['report_version']
    out_dict['Methylation_Classification_Final']=methyl_json['final_diagnosis']['methylation_class']
    out_dict['MGMT_Status']=methyl_json['final_diagnosis']['mgmt_status']

    levels = ["Superfamily", "Family", "Class", "Subclass"]

    if type is None:
        if 'IGM' in out_dict['Methylation_Version']:
            type = "IGM"
        elif "v12" in out_dict['Methylation_Version']:
            type = "v12"
        else:
            type = "v11"

    if type == "v11":
        print(f'This is a v11 methyl report JSON. Its version is {out_dict['Methylation_Version']}')
        # Detailed info for the v11 methylation are obtained from the rawdata JSON instead.
        pass
    elif type == "IGM":
        # This is an IGM classifier output, which has newer formatting.
        category = "Unclassified"
        level = "Unclassified"
        classifications = methyl_json['results']
        for i in classifications:
            category_string = i['category'].replace("Super Family","Superfamily")
            if category_string in levels:
                out_dict[f"Methylation_{category_string}"] = i['predictedClassification']
                out_dict[f"Methylation_{category_string}_Score"] = float(i["classifierScore"])
                if out_dict[f"Methylation_{category_string}_Score"] >= 0.80:
                    category = out_dict[f"Methylation_{category_string}"]
                    level = category_string
            else:
                print(f"Unknown methylation level: {i['category']} ({i['predictedClassification']})")
            
            out_dict["Methylation_Prediction_Category"] = category
            out_dict["Methylation_Prediction_Level"] = level
    elif type == "v12":

        methyl_scores = methyl_json['predicted_classification_classifier_scores']
        levels = ["Superfamily", "Family", "Class", "Subclass"]
        classifications = {}

        prev_level = None
        for i in methyl_scores:
            # Goes through the different category levels.
            # Tries to fix common parsing errors in the process
            category_string = i['category'].replace("Super Family","Superfamily")

            found = False
            for j in levels:
                if not found and category_string.startswith(j):
                    found = True
                    prev_level = j
                    category_string = category_string.replace(f"{j} ","")
                    score = i['score']
                    classifications[j]={"category":category_string, "score":score}
            if not found:
                if prev_level is not None:
                    for j in levels:
                        if j in category_string:
                            cat_string_parts = category_string.split(j)
                            category_string = cat_string_parts[0]
                            if j not in classifications:
                                classifications[j]={"category":cat_string_parts[1],"score":None}

                    classifications[prev_level]["category"] = (classifications[prev_level]["category"] + ", " + category_string).replace(",,",",")
                    if classifications[prev_level]['score'] is None and i['score'] is not None:
                        classifications[prev_level]['score'] = i['score']

        for i in classifications:
            if classifications[i]["score"] is None:
                for n in classifications[i]['category'].split(" "):
                    if not n.startswith("0"):
                        continue

                    rescued_score = None
                    try:
                        rescued_score = float(n)
                    except:
                        pass
                    if rescued_score is not None:
                        classifications[i]['score'] = rescued_score
                        classifications[i]['category'] = category_string.replace(f"{n} ","")
                        #print(f"Rescued a score: {category_string} [{score}]")
                        break

        category = "Unclassified"
        level = "Unclassified"
        for i in classifications:
            #print(f"{i} : {classifications[i]}")
            out_dict[f"Methylation_{i}"] = classifications[i]['category']
            out_dict[f"Methylation_{i}_Score"] = classifications[i]['score']
            if out_dict[f"Methylation_{i}_Score"] is None:
                print(f"{methyl_json['subject_id']} Methylation level {i} had missing score")
                out_dict[f"Methylation_{i}_Score"] = ""
            elif out_dict[f"Methylation_{i}_Score"] >= 0.90:
                category = out_dict[f"Methylation_{i}"]
                level = i
        out_dict["Methylation_Prediction_Category"] = category
        out_dict["Methylation_Prediction_Level"] = level
    else:
        print(f"Mode {type} is invalid. Skipping.")
        pass

    return out_dict

def parse_methyl_rawdata_json(methyl_data_json:dict, v11:bool=True, out_dict:dict={}):
    if methyl_data_json is None:
        return out_dict
    
    if v11:
        if "v11" not in methyl_ref:
            prep_methyl_ref_v11()

        print("Parsing a v11 raw data JSON.")
        family_id = ""
        family_score = ""
        class_id = ""
        class_score = ""

        if "family_data" in methyl_data_json:
            family_id = methyl_data_json["family_data"][0]["methylation_family"].strip().replace(",","").replace("/"," ").replace("  "," ").replace(" ","_").replace("MCF_","").replace("MTF_","").replace("MTGF_","")
            while "__" in family_id:
                family_id = family_id.replace("__","_")
            family_score = methyl_data_json["family_data"][0]["family_score"]
            if family_id in methyl_ref["v11"]:
                print(f"Family: {family_id} -> {methyl_ref["v11"][family_id]}")
                family_id = methyl_ref["v11"][family_id]
            else:
                print(f"Family {family_id} not found")

        if "class_data" in methyl_data_json:
            class_id = methyl_data_json["class_data"][0]["methylation_class"].strip().replace(",","").replace("/"," ").replace("  "," ").replace(" ","_")
            while "__" in class_id:
                class_id = class_id.replace("__","_")
            print(class_id)
            class_score = methyl_data_json["class_data"][0]["class_score"]
            if class_id in methyl_ref["v11"]:
                print(f"Class: {class_id} -> {methyl_ref['v11'][class_id]}")
                class_id = methyl_ref["v11"][class_id]
            else:
                print(f"Class {class_id} not found")

        out_dict["Methylation_Family"]=family_id
        out_dict["Methylation_Family_Score"]=family_score
        out_dict["Methylation_Class"]=class_id
        out_dict["Methylation_Class_Score"]=class_score

        #print(f"Family: {out_dict["Methylation_Family"]} [{out_dict["Methylation_Family_Score"]}] | Class: {out_dict["Methylation_Class"]} [{out_dict["Methylation_Class_Score"]}]")

        if "MGMT_Status" not in out_dict:
            out_dict["MGMT_Status"]=methyl_data_json["mgmt_methylation_data"][0]["mgmt_methylation_status"]
        
    else:
        # Should not need to parse v12 raw data files.
        pass

    return out_dict


### ARCHER Fusion

def parse_archer_json(archer_json, out_dict:dict={}):
    # IGM Archer gene fusion & intragenic break detection
    if archer_json is None:
        return out_dict
    out_dict = parse_molecular_generic(archer_json, out_dict)

    if "DKFZ" in archer_json['report_version']:
        print("Something has gone wrong. Why is this methylation file being processed as Archer?")

    out_dict["Archer_Version"] = archer_json['report_version']
    out_dict['Archer_Result_Tier1-2'] = 'Negative'
    out_dict['Archer_Result_Tier3'] = 'Negative'

    if 'fusion_tier_one_or_two_result' in archer_json:
        results = []
        for i in archer_json['fusion_tier_one_or_two_result']['variants']:
            results.append(i['gene_fusion'])
        out_dict['Archer_Tier1-2_Fusions']=";".join(results)
        if len(out_dict['Archer_Tier1-2_Fusions']) > 0:
            out_dict['Archer_Result_Tier1-2'] = 'Positive'
        if 'summary' in archer_json['fusion_tier_one_or_two_result']:
            out_dict['Archer_Blurb_Tier1-2'] = " ".join(archer_json['fusion_tier_one_or_two_result']['summary'])
    if 'single_tier_one_or_two_result' in archer_json:
        results = []
        for i in archer_json['single_tier_one_or_two_result']['variants']:
            results.append(i['breakpoint1']['gene'])
        out_dict['Archer_Tier1-2_Intragenic']=";".join(results)
        if len(out_dict['Archer_Tier1-2_Intragenic']) > 0:
            out_dict['Archer_Result_Tier1-2'] = 'Positive'
        if 'summary' in archer_json['single_tier_one_or_two_result']:
            out_dict['Archer_Blurb_Tier1-2_Intragenic'] = " ".join(archer_json['single_tier_one_or_two_result']['summary'])
    if 'fusion_tier_three_result' in archer_json:
        results = []
        for i in archer_json['fusion_tier_three_result']['variants']:
            results.append(i['gene_fusion'])
        out_dict['Archer_Tier3_Fusions']=";".join(results)
        if len(out_dict['Archer_Tier3_Fusions']) > 0:
            out_dict['Archer_Result_Tier3'] = 'Positive'
        if 'summary' in archer_json['fusion_tier_three_result']:
            out_dict['Archer_Blurb_Tier3'] = " ".join(archer_json['fusion_tier_three_result']['summary'])
    if 'single_tier_three_result' in archer_json:
        results = []
        for i in archer_json['single_tier_three_result']['variants']:
            results.append(i['breakpoint1']['gene'])
        out_dict['Archer_Tier3_Intragenic']=";".join(results)
        if len(out_dict['Archer_Tier3_Intragenic']) > 0:
            out_dict['Archer_Result_Tier3'] = 'Positive'
        if 'summary' in archer_json['single_tier_three_result']:
            out_dict['Archer_Blurb_Tier3_Intragenic'] = " ".join(archer_json['single_tier_three_result']['summary'])

    return out_dict

# Passes samples through the above parsers in sequence

def parse_sample_jsons(sample_jsons:dict):
    out_dict = {}

    if "cog" in sample_jsons:
        out_dict = parse_cog_json(sample_jsons["cog"], out_dict)
    if "tumor_normal" in sample_jsons:
        out_dict = parse_tumor_normal_json(sample_jsons["tumor_normal"], out_dict)
    if ("methyl_igm" in sample_jsons and sample_jsons["methyl_igm"] is not None) or ("methyl_v12" in sample_jsons and sample_jsons["methyl_v12"] is not None):
        if "methyl_igm" in sample_jsons and sample_jsons["methyl_igm"] is not None:
            out_dict = parse_methyl_json(sample_jsons["methyl_igm"], "IGM", out_dict)
        elif "methyl_v12" in sample_jsons and sample_jsons["methyl_v12"] is not None:
            out_dict = parse_methyl_json(sample_jsons["methyl_v12"], "v12", out_dict)
    elif ("methyl_v11" in sample_jsons and sample_jsons["methyl_v11"] is not None) or ("methyl_v11_raw" in sample_jsons and sample_jsons["methyl_v11_raw"] is not None):
        if "methyl_v11_raw" in sample_jsons and sample_jsons["methyl_v11_raw"] is not None:
            out_dict = parse_methyl_rawdata_json(sample_jsons["methyl_v11_raw"], True, out_dict)
        if "methyl_v11" in sample_jsons and sample_jsons["methyl_v11"] is not None:
            out_dict = parse_methyl_json(sample_jsons["methyl_v11"], "v11", out_dict)
    if "archer_fusion" in sample_jsons:
        out_dict = parse_archer_json(sample_jsons["archer_fusion"], out_dict)
    #print(out_dict)
    return out_dict

# Main function

def run_json_parser(args: argparse.Namespace) -> None:
    import pandas as pd
    import os

    json_dirs = args.input_json_dirs.split(",") #"/sbgenomics/project-files/*/" #"/Users/glw001/Projects/MCI_Report_JSONs/MCI_4-5-2024/"
    out_dir = os.path.dirname(args.output_prefix)
    #file_prefix = os.path.basename(args.output_prefix)
    excel_out = os.path.join(out_dir, f"{args.output_prefix}.xlsx")
    json_out = os.path.join(out_dir, f"{args.output_prefix}.json")

    debug=True

    blank_field_placeholder = args.blank_field_indicator
    mci_dict_reference = args.data_dict_reference

    if debug:
        print("Processing data dictionary...")
    data_dict_table = pd.read_csv(mci_dict_reference, delimiter="\t", header=0,index_col = 0,keep_default_na=False)
    if debug:
        print(data_dict_table)

    if debug:
        print("Reading in data JSONs...")
    json_dicts = sort_jsons(json_dirs, debug)

    if debug:
        print(f"Got {len(json_dicts)} samples' data.")

    data = []

    for i in json_dicts:
        #print(i)
        sample_jsons = json_dicts[i]
        out_dict = replace_blank_fields(parse_sample_jsons(sample_jsons), blank_field_placeholder)
        out_dict['Sample']=i
        data.append(out_dict)

    data = standardize_variant_notation(data, debug=debug)
    
    data = standardize_methylation_class(data, debug=debug)

    out_df = pd.DataFrame(data, columns=list(data_dict_table['Term']))

    if args.output_type.lower() in ['excel','both']:
        with pd.ExcelWriter(excel_out) as writer:  
            out_df.to_excel(writer, sheet_name='MCI JSON Data',index=False)
            data_dict_table.to_excel(writer, sheet_name='Data Dictionary',index=False)
            print(f"Excel sheet writen to: {excel_out}")
    
    if args.output_type.lower() in ['json','both']:
        json_formatted_data = {}
        for i in data:
            json_formatted_data[i['Sample']]=i
        datadict_as_dict = data_dict_table.to_dict()
        json_formatted_datadict = {}
        for i in range(len(datadict_as_dict['Term'])):
            term = datadict_as_dict['Term'][i]
            definition = datadict_as_dict['Definition'][i]
            source = datadict_as_dict['JSON Source'][i]
            note = datadict_as_dict['Notes'][i]
            rave_id = datadict_as_dict['RAVE Identifier / JSON Field'][i]
            json_formatted_datadict[term] = {"Definition":definition, "Source":source, "Note":note, "RAVE Identifier or JSON Field":rave_id}

        json_dict = {'data':json_formatted_data, 'dictionary':json_formatted_datadict}
        with open(json_out,'w') as json_file:
            json_file.write(json.dumps(json_dict))
        print(f"JSON written to: {json_out}")

def mci_json_argparser():
    parser = argparse.ArgumentParser(
        prog=mci_json_argparser.__name__, description="")
    parser.add_argument(
        '--input-json-dirs', required=True,
        help="String of comma separated directory paths where JSON files are located.")
    parser.add_argument(
        '--output-type',
        help="Format of output data file(s).", choices=["Excel","JSON","Both"], default="Both")
    parser.add_argument(
        '--data-dict-reference', required=False, default="./mci_data_dict.txt",
        help="Path to text file containing data dictionary to include in Excel outputs.")
    parser.add_argument(
        '--blank-field-indicator', required=False, default=".",
        help="Fields that are present but blank (as opposed to missing) will be indicated with the specified string.")
    parser.add_argument(
        '--output-prefix', type=str, required=True,
        help="Processed data output prefix.")

    return parser

def main():
    parser = mci_json_argparser()
    args = parser.parse_args()
    run_json_parser(args)

if __name__ == "__main__":
    main()