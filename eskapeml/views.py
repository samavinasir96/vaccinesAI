from pathlib import Path
from django.shortcuts import render
import subprocess
import time
import zipfile
import iedb
import joblib
import numpy as np
import pandas as pd
import os
import sys    
from django.http import HttpRequest, HttpResponse, HttpResponseRedirect, JsonResponse
from django.views.decorators.csrf import csrf_exempt
from django.shortcuts import redirect, render
import pandas as pd
import os
import sys
import logging
import zipfile
from Bio import SeqIO
import biolib
from tensorflow.keras.models import load_model
from django.core.cache import cache

logger = logging.getLogger(__name__)

# Configure logging
def setup_logger():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    log_file_path = os.path.join(script_dir, "vacsol_ml.log")
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file_path),
            logging.StreamHandler()
        ]
    )
setup_logger()

def get_progress(request):
    """Endpoint to get both progress percentage and logs"""
    progress = cache.get(f'progress_{request.session.session_key}', 0)
    logs = cache.get(f'logs_{request.session.session_key}', [])
    print(f"Progress: {progress}, Logs: {logs}")
    return JsonResponse({
        'progress': progress,
        'logs': logs
    })

def progress_callback(request, message, progress):
    """Callback to update progress and store logs"""

    print(f"Setting progress: {progress}% - {message}")

    cache.set(f'progress_{request.session.session_key}', progress, timeout=3600)    
    logs = cache.get(f'logs_{request.session.session_key}', [])
    logs.append(f"{progress}% - {message}")
    cache.set(f'logs_{request.session.session_key}', logs[-10:], timeout=3600)  # Keep last 10 messages
    
    # Also log to file
    logger.info(f"{progress}% - {message}")
    

@csrf_exempt
def home(request):
    return render (request, 'index.html')

def faqs(request):
    return render (request, 'faqs.html')

def glossary(request):
    return render (request, 'glossary.html')

def get_latest_logs(request):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    log_file_path = os.path.join(script_dir, "vacsol_ml.log")
    with open(log_file_path, 'r') as log_file:
        latest_logs = log_file.readlines()
    return JsonResponse({'logs': latest_logs})

def upload_sequence(request):
    if request.method == "POST":
        if not request.session.session_key:
            request.session.create()

        try:
            cache.set(f'progress_{request.session.session_key}', 0, timeout=3600)
            cache.set(f'logs_{request.session.session_key}', [], timeout=3600)

            sequence = request.POST.get("sequence")
            #sequence = sequence.replace("\n", "")
            file = request.FILES.get("file")

            if sequence:
                script_dir = os.path.dirname(os.path.abspath(__file__))
                file_path = os.path.join(script_dir, "sequences.fasta")

                with open(file_path, "w") as f:
                    f.write(sequence)
            elif file:
                script_dir = os.path.dirname(os.path.abspath(__file__))
                file_path = os.path.join(script_dir, "sequences.fasta")
                with open(file_path, "wb") as f:
                    for chunk in file.chunks():
                        f.write(chunk)
            
            else:
                cache.delete(f'progress_{request.session.session_key}')
                return JsonResponse(
                    {"status": "error", "message": "No sequence or file provided"},
                    status=400
                )
            
            cache.set(f'progress_{request.session.session_key}', 0, timeout=3600)
            cache.set(f'logs_{request.session.session_key}', ["Starting analysis..."], timeout=3600)
            
            calculate_features(request, file_path, lambda message, progress: progress_callback(request, message, progress))

            return JsonResponse({"status": "processing", "redirect": "/eskapeml_results"})
        
        except Exception as e:
            logger.error(f"Upload error: {str(e)}")
            cache.delete(f'progress_{request.session.session_key}')
            return JsonResponse({"status": "error", "message": str(e)}, status=500)

        # return redirect("eskapeml_results")
    return render(request, 'upload_vacsolml.html')

def processing_view(request):
    return render(request, 'processing.html')   

def calculate_features(request, file_path, progress_callback):
    
    total_steps = 20
    completed_steps = 0

    script_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(script_dir, "sequences.fasta")
    file_path = os.path.abspath(file_path)
    MERGED_DATAFRAMES = []
    mhci_prob_scores = []
    mhci_ranks = []
    mhcii_ranks = []
    mhcii_scores = []
    b_cells_prob_scores = []
    surface_probs = []
    antigenicities = []
    header_list = []

    if progress_callback:
        completed_steps += 1
        progress_callback(message="Sequence(s) uploaded successfully", progress=int((completed_steps/total_steps)*100))

    try:
        if progress_callback:
            completed_steps += 1
            progress_callback(message="Physiochemical parameters analysis in progress", progress=int((completed_steps/total_steps)*100))

        # Get base directory (adjust if your structure differs)
        BASE_DIR = Path(__file__).resolve().parent.parent

        # Define absolute paths
        IFEATURE_PATH = BASE_DIR / "Data" / "iFeature" / "iFeature.py"
        INPUT_FILE = BASE_DIR / "eskapeml" / "sequences.fasta"
        OUTPUT_DIR = BASE_DIR / "eskapeml" / "Temp_Results"

        # Create output directory if it doesn't exist
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

        commands = [
            ["python", str(IFEATURE_PATH), "--file", str(INPUT_FILE), "--type", "Moran", "--out", str(OUTPUT_DIR/"ifeature3.tsv")],
            ["python", str(IFEATURE_PATH), "--file", str(INPUT_FILE), "--type", "Geary", "--out", str(OUTPUT_DIR/"ifeature4.tsv")],
            ["python", str(IFEATURE_PATH), "--file", str(INPUT_FILE), "--type", "NMBroto", "--out", str(OUTPUT_DIR/"ifeature5.tsv")],
            ["python", str(IFEATURE_PATH), "--file", str(INPUT_FILE), "--type", "CTriad", "--out", str(OUTPUT_DIR/"ifeature9.tsv")],
            ["python", str(IFEATURE_PATH), "--file", str(INPUT_FILE), "--type", "KSCTriad", "--out", str(OUTPUT_DIR/"ifeature10.tsv")],
            ["python", str(IFEATURE_PATH), "--file", str(INPUT_FILE), "--type", "SOCNumber", "--out", str(OUTPUT_DIR/"ifeature11.tsv")],
            ["python", str(IFEATURE_PATH), "--file", str(INPUT_FILE), "--type", "QSOrder", "--out", str(OUTPUT_DIR/"ifeature12.tsv")],
            ["python", str(IFEATURE_PATH), "--file", str(INPUT_FILE), "--type", "CKSAAP", "--out", str(OUTPUT_DIR/"ifeature15.tsv")],
            ["python", str(IFEATURE_PATH), "--file", str(INPUT_FILE), "--type", "CKSAAGP", "--out", str(OUTPUT_DIR/"ifeature19.tsv")],
        ]

        for cmd in commands:
            try:
                result = subprocess.run(
                    cmd,
                    check=True,
                    capture_output=True,
                    text=True
                )
                print(f"Success: {' '.join(cmd)}")
                print(result.stdout)
            except subprocess.CalledProcessError as e:
                print(f"Error running: {' '.join(cmd)}")
                raise
        

        subprocess.run('python "Data/iFeature/iFeature.py" --file eskapeml/sequences.fasta --type Moran --out eskapeml/Temp_Results/ifeature3.tsv')
        subprocess.run('python "Data/iFeature/iFeature.py" --file eskapeml/sequences.fasta --type Geary --out eskapeml/Temp_Results/ifeature4.tsv')
        subprocess.run('python "Data/iFeature/iFeature.py" --file eskapeml/sequences.fasta --type NMBroto --out eskapeml/Temp_Results/ifeature5.tsv')
        subprocess.run('python "Data/iFeature/iFeature.py" --file eskapeml/sequences.fasta --type CTriad --out eskapeml/Temp_Results/ifeature9.tsv')
        subprocess.run('python "Data/iFeature/iFeature.py" --file eskapeml/sequences.fasta --type KSCTriad --out eskapeml/Temp_Results/ifeature10.tsv')
        subprocess.run('python "Data/iFeature/iFeature.py" --file eskapeml/sequences.fasta --type SOCNumber --out eskapeml/Temp_Results/ifeature11.tsv')
        subprocess.run('python "Data/iFeature/iFeature.py" --file eskapeml/sequences.fasta --type QSOrder --out eskapeml/Temp_Results/ifeature12.tsv')
        subprocess.run('python "Data/iFeature/iFeature.py" --file eskapeml/sequences.fasta --type CKSAAP --out eskapeml/Temp_Results/ifeature15.tsv')
        subprocess.run('python "Data/iFeature/iFeature.py" --file eskapeml/sequences.fasta --type CKSAAGP --out eskapeml/Temp_Results/ifeature19.tsv')
        
        if progress_callback:
            completed_steps += 1
            progress_callback(message="Physicochemical analysis completed", progress=int((completed_steps/total_steps)*100))

        df3 = pd.read_table ('eskapeml/Temp_Results/ifeature3.tsv')
        df4 = pd.read_table('eskapeml/Temp_Results/ifeature4.tsv')
        df5 = pd.read_table('eskapeml/Temp_Results/ifeature5.tsv')
        df9 = pd.read_table ('eskapeml/Temp_Results/ifeature9.tsv')
        df10 = pd.read_table('eskapeml/Temp_Results/ifeature10.tsv')
        df11 = pd.read_table('eskapeml/Temp_Results/ifeature11.tsv')
        df12 = pd.read_table ('eskapeml/Temp_Results/ifeature12.tsv')
        df15 = pd.read_table('eskapeml/Temp_Results/ifeature15.tsv')
        df19 = pd.read_table ('eskapeml/Temp_Results/ifeature19.tsv')

        MERGED_DATAFRAMES = []
        df_all = [df3, df4, df5, df9, df10, df11, df12, df15, df19]
        con1 = pd.concat(df_all, axis="columns")
        con1 = con1.loc[:, ~con1.columns.duplicated()] 
        MERGED_DATAFRAMES.append(con1)

        if progress_callback:
            completed_steps += 1
            progress_callback(message="Creating Dataset", progress=int((completed_steps/total_steps)*100))

        con1.to_csv ('eskapeml/Temp_Results/eskape_ifeature.csv', index=False)

        df_final = pd.read_csv('eskapeml/Temp_Results/eskape_ifeature.csv')
        #columns_to_keep = ['#','g3.g4.g3', 'AN.gap0', 'EM.gap0', 'LV.gap0', 'SV.gap0', 'QG.gap1', 'TI.gap1', 'FM.gap2', 'GQ.gap2', 'WC.gap2', 'RY.gap3', 'GD.gap4', 'YG.gap4', 'GS.gap5', 'PF.gap5', 'TR.gap5', 'TV.gap5', 'VI.gap5', 'ALF', 'AMP', 'APR', 'AVT', 'EFR', 'ELA', 'EYA', 'FKR', 'GGQ', 'GLR', 'GVP', 'HKG', 'IDD', 'IDR', 'IYQ', 'KLQ', 'KLS', 'LEK', 'NLT', 'NRY', 'NSG', 'PPV', 'RSL', 'RSS', 'SGE', 'SGH', 'TSV', 'VAG']
        #final_df = df_final[columns_to_keep]
        final_df = df_final.rename(columns={'#': 'Protein_ID'})
        final_df.to_csv('eskapeml/Temp_Results/ifeatures.csv', index=False)

        #features_file for users
        script_dir = sys.path[0]
        results_folder_path = os.path.join(script_dir, "eskapeml/Analysis_Results")

        #physicochemical features
        physico_filename = 'physiochemical_parameters.csv'
        physico_result_path = os.path.join(results_folder_path, physico_filename)
        final_df.to_csv(physico_result_path, index=False)   

        fasta_file_path = "eskapeml/sequences.fasta"

        # Read the FASTA file
        sequences = SeqIO.to_dict(SeqIO.parse(fasta_file_path, "fasta"))

        # Loop through each sequence and send POST requests
        for header, sequence_record in sequences.items():
            sequence = str(sequence_record.seq)
            header_list.append(header)

            if progress_callback:
                completed_steps += 1
                progress_callback(message="MHC class I epitope prediction in progress", progress=int((completed_steps/total_steps)*100))

            # Send POST request to MHC class I peptide binding prediction tool:
            mhci_res1 = iedb.query_mhci_binding(method="recommended", sequence= sequence, allele="HLA-A*02:01", length="9")
            mhci_res2 = iedb.query_mhci_binding(method="recommended", sequence= sequence, allele="HLA-A*01:02", length="9")

            #concat mhci-res
            mhci_res = pd.concat([mhci_res1, mhci_res2], axis=0)
            time.sleep(3)

        for header, sequence_record in sequences.items():
            sequence = str(sequence_record.seq)
            header_list.append(header)

            if progress_callback:
                completed_steps += 1
                progress_callback(message="MHC class II epitope prediction in progress", progress=int((completed_steps/total_steps)*100))

            # Send POST request to MHC class II peptide binding prediction tool:
            mhcii_res1 = iedb.query_mhcii_binding(method="recommended", sequence= sequence, allele="HLA-DRB1*01:01", length=None)
            mhcii_res2 = iedb.query_mhcii_binding(method="recommended", sequence= sequence, allele="HLA-DRB1*01:02", length=None)
            
            mhcii_res = pd.concat([mhcii_res1, mhcii_res2], axis=0)
            time.sleep(3)

        for header, sequence_record in sequences.items():
            sequence = str(sequence_record.seq)
            header_list.append(header)
            
            if progress_callback:
                completed_steps += 1
                progress_callback(message="B-cell epitope prediction in progress", progress=int((completed_steps/total_steps)*100))

            bcell_res = iedb.query_bcell_epitope(method="Bepipred", sequence= sequence, window_size=9)

            if progress_callback:
                completed_steps += 1
                progress_callback(message="Epitope Analysis Completed", progress=int((completed_steps/total_steps)*100))

            time.sleep(3)

        for header, sequence_record in sequences.items():
            sequence = str(sequence_record.seq)
            header_list.append(header)
            
            if progress_callback:
                completed_steps += 1
                progress_callback(message="Adhesion Probability prediction in progress", progress=int((completed_steps/total_steps)*100))

            sprob_res = iedb.query_bcell_epitope(method="Emini", sequence= sequence, window_size=9)

            if progress_callback:
                completed_steps += 1
                progress_callback(message="Adhesion Probability Analysis Completed", progress=int((completed_steps/total_steps)*100))

            time.sleep(3)

        for header, sequence_record in sequences.items():
            sequence = str(sequence_record.seq)
            header_list.append(header)

            if progress_callback:
                completed_steps += 1
                progress_callback(message="Antigenicity prediction in progress", progress=int((completed_steps/total_steps)*100))

            antigenicity_1 = iedb.query_bcell_epitope(method="Kolaskar-Tongaonkar", sequence= sequence, window_size=9)

            if progress_callback:
                completed_steps += 1
                progress_callback(message="Antigenicity Prediction Completed", progress=int((completed_steps/total_steps)*100))

            time.sleep(3)

            # Getting means - mhci
            mhci_res['score'] = pd.to_numeric(mhci_res['score'], errors='coerce')
            mhci_res['percentile_rank'] = pd.to_numeric(mhci_res['percentile_rank'], errors='coerce')
                
            mhci_res_filtered = mhci_res[mhci_res["percentile_rank"] <= 0.5]
                
            if not mhci_res_filtered.empty:
                df_score_mhci = mhci_res_filtered["score"].mean()
                df_rank_mhci = mhci_res_filtered["percentile_rank"].mean()
                mhci_epitopes = mhci_res_filtered['peptide'].values
            else:
                df_score_mhci= 0
                df_rank_mhci= 0 

            # Getting means - mhcii
            mhcii_res['score'] = pd.to_numeric(mhcii_res['score'], errors='coerce')
            mhcii_res['rank'] = pd.to_numeric(mhcii_res['rank'], errors='coerce')
            
            mhcii_res_filtered = mhcii_res[mhcii_res["rank"] <= 1]

            if not mhcii_res_filtered.empty:
                rank_mhcii = mhcii_res_filtered["rank"].mean()
                score_mhcii = mhcii_res_filtered["score"].mean()
                mhcii_epitopes = mhcii_res_filtered['peptide'].values      
            else:
                rank_mhcii= 0
                score_mhcii= 0 

            # Getting means - bcells
            bcell_res['Score'] = pd.to_numeric(bcell_res['Score'], errors='coerce')
            df_bcells_final = bcell_res["Score"].mean()

            # Getting means - s_probability
            sprob_res['Score'] = pd.to_numeric(sprob_res['Score'], errors='coerce')
            df_sprob_final = sprob_res["Score"].mean()

            # Getting means - antigenicity(scale1)
            antigenicity_1['Score'] = pd.to_numeric(antigenicity_1['Score'], errors='coerce')
            antigenicity_final = antigenicity_1["Score"].mean()

            if progress_callback:
                completed_steps += 1
                progress_callback(message="Updating Dataset", progress=int((completed_steps/total_steps)*100))

            add_scores2 = pd.read_csv('eskapeml/Temp_Results/ifeatures.csv')
            mhci_prob_scores.append(df_score_mhci)
            mhci_ranks.append(df_rank_mhci)
            mhcii_ranks.append(rank_mhcii)
            mhcii_scores.append(score_mhcii)
            b_cells_prob_scores.append(df_bcells_final)
            surface_probs.append(df_sprob_final)
            antigenicities.append(antigenicity_final)

            # Create a new DataFrame using the lists
            result_df = pd.DataFrame({
                #'Protein_ID': protein_ids,
                'antigenicity_1': antigenicities,
                'b_cells_probability_score': b_cells_prob_scores,
                'mhci_probability_score': mhci_prob_scores,
                'mhci_rank': mhci_ranks,
                'mhcii_rank': mhcii_ranks,
                'mhcii_score': mhcii_scores,
                'surface_probability': surface_probs,
            })
            result_df = result_df.round(6)
            merged_df = pd.concat([add_scores2, result_df], axis=1)
            merged_df.to_csv('eskapeml/Temp_Results/ifeatures_updated.csv', index=False)

            if progress_callback:
                completed_steps += 1
                progress_callback(message="Dataset Updated", progress=int((completed_steps/total_steps)*100))

        if progress_callback:
            completed_steps += 1
            progress_callback(message="Signal peptides prediction in progress", progress=int((completed_steps/total_steps)*100))

        signalp_6 = biolib.load('DTU/SignalP_6')
        signalp6_job = signalp_6.cli(args='--fastafile eskapeml/sequences.fasta --output_dir output')
        signalp6_job.save_files('eskapeml/Temp_Results/signalP')

        script_dir = sys.path[0]
        signalP_folder_path = os.path.join(script_dir, "eskapeml/Temp_Results/signalP")

        if progress_callback:
            completed_steps += 1
            progress_callback(message="Signal Peptide Analysis Completed", progress=int((completed_steps/total_steps)*100))

        #read signalp results:    
        sp_table_path = os.path.join(signalP_folder_path, "prediction_results.txt")
        df = pd.read_table(sp_table_path, sep="\t", header=None, skiprows=1)
        df.columns = ["ID", "Prediction", "OTHER", "SP(Sec/SPI)", "LIPO(Sec/SPII)", "TAT(Tat/SPI)", "TATLIPO(Sec/SPII)", "PILIN(Sec/SPIII)", "CS Position"]
        df.drop(index=df.index[0], axis=1, inplace=True)

        #extract the scores
        SP = df['SP(Sec/SPI)'].values
        LIPO = df['LIPO(Sec/SPII)'].values
        TAT = df['TAT(Tat/SPI)'].values
        TATLIPO = df['TATLIPO(Sec/SPII)'].values
        PILIN = df['PILIN(Sec/SPIII)'].values
        OTHER = df['OTHER'].values
        
        #features_file for users
        script_dir = sys.path[0]
        results_folder_path = os.path.join(script_dir, "eskapeml/Analysis_Results")

        #Signal peptides features
        threshold_signalp = [
                    'SignalP result threshold: If a value of any column is 1 or close to 1, the protein belongs to that category',
                ]
        threshold_signalp_df = pd.DataFrame({'Thresholds': threshold_signalp})

        signalp_data = {
        'signal_peptide_SP': SP,
        'signal_peptide_LIPO': LIPO,
        'signal_peptide_TATLIPO': TATLIPO,
        'signal_peptide_TAT': TAT,
        'signal_peptide_PILIN': PILIN,
        'signal_peptide_OTHERS': OTHER  
        }

        df_signalp = pd.DataFrame(signalp_data)

        biological_filename = 'biological_properties.csv'

        biological_df =  pd.DataFrame({
            'Protein_ID': header,
            'antigenicity_1': antigenicities,
            'b_cells_probability_score': b_cells_prob_scores,
            'mhci_probability_score': mhci_prob_scores,
            'mhci_rank': mhci_ranks,
            'mhcii_rank': mhcii_ranks,
            'mhcii_score': mhcii_scores,
            'surface_probability': surface_probs,
        })

        biological_df = biological_df.round(6)


        threshold_adhesion = [
            'Threshold for adhesion probability = 1 Or greater than 1',
        ]
        threshold_adesion_df = pd.DataFrame({'Thresholds': threshold_adhesion})

        
        df_mhci_epitopes = pd.DataFrame({
                                            'Epitope_Type': ['mhci'] * len(mhci_epitopes),
                                            'Epitope_Sequence': mhci_epitopes})
        
        df_mhcii_epitopes = pd.DataFrame({
                                            'Epitope_Type': ['mhcii'] * len(mhcii_epitopes),
                                            'Epitope_Sequence': mhcii_epitopes})
        
        comments = [
            'Threshold for strong "MHC Class I" binder: % Rank = Lower than 0.5',
            'Threshold for strong "MHC Class II" binder: % Rank = Lower than 1',
            'Threshold for strong "B cell epitope": Score = Greater than 0.5',
            'Antigenicity = Value close to 1 or Greater than 1 indicates high antigenicity'
        ]
        
        biological_result_path = os.path.join(results_folder_path, biological_filename)
        
        comments_df = pd.DataFrame({'Thresholds': comments})

        df_biological = pd.concat([biological_df, df_signalp, df_mhci_epitopes, df_mhcii_epitopes, comments_df, threshold_signalp_df], ignore_index=True)
            
        df_biological.to_csv(biological_result_path, index=False)

        if progress_callback:
            completed_steps += 1
            progress_callback(message="Dataset Updated", progress=int((completed_steps/total_steps)*100))

        add_scores = pd.read_csv('eskapeml/Temp_Results/ifeatures_updated.csv')
        add_scores["signal_peptide_SP"] = SP
        add_scores["signal_peptide_LIPO"] = LIPO
        add_scores["signal_peptide_TAT"] = TAT
        add_scores["signal_peptide_TATLIPO"] = TATLIPO
        add_scores["signal_peptide_PILIN"] = PILIN
        add_scores["signal_peptide_OTHER"] = OTHER

        add_scores.to_csv('eskapeml/Temp_Results/finalfeatures_updated.csv', index=False)

        final = pd.read_csv('eskapeml/Temp_Results/finalfeatures_updated.csv')

        Protein_ID = final['Protein_ID']
        finalDF = final.drop(labels=["Protein_ID"], axis=1)

        if progress_callback:
            completed_steps += 1
            progress_callback(message="Encoding Selected ESKAPE Pathogen", progress=int((completed_steps/total_steps)*100))

        pathogen_encoding = {
            'Staphylococcus aureus': 5,
            'Enterococcus faecium': 2,
            'Pseudomonas aeruginosa': 4,
            'Acinetobacter baumannii': 0,
            'Klebsiella pneumoniae': 3,
            'Enterobacter': 1
        }

        if request.method == 'POST':
            selected_pathogen = request.POST.get('pathogen')
            
            encoded_pathogen = pathogen_encoding.get(selected_pathogen, -1)  # Default to -1 if not found
            
            if encoded_pathogen != -1:

                finalDF["Organism.1"] = encoded_pathogen

                finalDF.to_csv('eskapeml/Temp_Results/finalfeatures_pathogen.csv', index=False)
            else:
                print("Error: Selected pathogen not found in the encoding dictionary.")


        scaler_path = "eskapeml/eskape_scaler.pkl"

        if progress_callback:
            completed_steps += 1
            progress_callback(message="Scaling Dataset", progress=int((completed_steps/total_steps)*100))

        dataset = pd.read_csv('eskapeml/Temp_Results/finalfeatures_pathogen.csv')

        with open(scaler_path, "rb") as scaler:
            scale_dataset = joblib.load(scaler)    
            rescaled=scale_dataset.transform(dataset)
        
        rescaleData = pd.DataFrame(rescaled, columns=finalDF.columns)
        #rescaleData.insert(0, 'Protein_ID', Protein_ID)
        rescaleData.to_csv('eskapeml/Temp_Results/final_dataset.csv', index=False)

        if progress_callback:
            completed_steps += 1
            progress_callback(message="Features Annotation Completed - Saving Results", progress=int((completed_steps/total_steps)*100))
    
    except Exception as e:
        if progress_callback:
            progress_callback(f"Error: {str(e)}", 0)
        raise
    

def update_progress(request, message, step):
    """Update progress and log message"""

    progress = int((step / 20) * 100)
    cache.set(f'progress_{request.session.session_key}', progress, timeout=3600)
    
    logger = logging.getLogger(__name__)
    logger.info(message, extra={'progress': progress})
    
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    log_file_path = os.path.join(script_dir, "vacsol_ml.log")
    file_handler = logging.FileHandler(log_file_path)
    file_handler.setLevel(logging.INFO)
    logger.addHandler(file_handler)

    with open(log_file_path, 'r') as log_file:
        logs = log_file.readlines()

    # Clear the log file after reading the messages
    with open(log_file_path, 'w') as log_file:
        log_file.truncate(0)

    return render(request, 'upload_vacsolml.html', {'logs': logs})

    
import csv

def get_results(request):
    final = pd.read_csv('eskapeml/Temp_Results/finalfeatures_updated.csv')

    try:
        script_dir = sys.path[0]

        feature_cols = pd.read_csv('eskapeml/Temp_Results/final_dataset.csv')

        encoder = load_model('eskapeml/encoder_model.keras')
        encoded_data = encoder.predict(feature_cols)

        model_path = os.path.join(script_dir, "eskapeml/eskape_ensemble_model.pkl")

        with open(model_path, "rb") as model:
            ensembl_model = joblib.load(model)

        rf_model = ensembl_model['Random_Forest']
        gb_model = ensembl_model['Gradient_Boosting']
        lr_model = ensembl_model['Logistic_Regression']

        rf_predictions = rf_model.predict(encoded_data)
        gb_predictions = gb_model.predict(encoded_data)
        lr_predictions = lr_model.predict(encoded_data)

        rf_proba = rf_model.predict_proba(encoded_data)
        gb_proba = gb_model.predict_proba(encoded_data)
        lr_proba = lr_model.predict_proba(encoded_data)

        combined_predictions = np.where((rf_predictions + gb_predictions + lr_predictions) >= 1, 1, 0)
        combined_proba = (rf_proba + gb_proba + lr_proba) / 3
        print(combined_proba)


        proba_class_1 = combined_proba[:, 0]  
        proba_class_0 = combined_proba[:, 1]
        print(proba_class_1)


        csv_file = 'eskapeml/Temp_Results/finalfeatures_pathogen.csv'

        features_view = pd.read_csv(csv_file)

        results = pd.DataFrame({"Protein_ID": final['Protein_ID'], "Prediction": combined_predictions, "Probability_Class_1": proba_class_1, "Probability_Class_0": proba_class_0,
                                "antigenicity_1": features_view["antigenicity_1"], "b_cells_probability_score": features_view["b_cells_probability_score"],
                                "mhci_probability_score": features_view["mhci_probability_score"], "mhci_rank": features_view["mhci_rank"], "mhcii_rank": features_view["mhcii_rank"], "mhcii_score": features_view["mhcii_score"],
                                 "surface_probability": features_view["surface_probability"], "signal_peptide_SP": features_view["signal_peptide_SP"], "signal_peptide_LIPO": features_view["signal_peptide_LIPO"], "signal_peptide_TAT": features_view["signal_peptide_TAT"], "signal_peptide_TATLIPO": features_view["signal_peptide_TATLIPO"], "signal_peptide_PILIN": features_view["signal_peptide_PILIN"],
                                })

        results_data = [
            {"Protein_ID": row["Protein_ID"], "Prediction": int(row["Prediction"]), "Probability_Class_1": row["Probability_Class_1"], "Probability_Class_0": row["Probability_Class_0"],
             "antigenicity_1": row["antigenicity_1"], "b_cells_probability_score": row["b_cells_probability_score"],
             "mhci_probability_score": row["mhci_probability_score"], "mhci_rank": row["mhci_rank"], "mhcii_rank": row["mhcii_rank"], "mhcii_score": row["mhcii_score"],
             "surface_probability": row["surface_probability"], "signal_peptide_SP": row["signal_peptide_SP"], "signal_peptide_LIPO": row["signal_peptide_LIPO"], "signal_peptide_TAT": row["signal_peptide_TAT"], "signal_peptide_TATLIPO": row["signal_peptide_TATLIPO"], "signal_peptide_PILIN": row["signal_peptide_PILIN"],
             }
            for _, row in results.iterrows()
        ]

        return render(request, "results_vacsolml.html", {"eskapeml_results": results_data})

    except Exception as e:
        print(f"Error analyzing sequence: {e}")
        return JsonResponse({"status": "Error analyzing sequence"}, status=500)

def display_features(request):
    selected_columns = ['antigenicity_1','b_cells_probability_score','mhci_probability_score','mhci_rank,mhcii_rank']  # Specify the columns you want to display
    csv_file = 'eskapeml/Temp_Results/finalfeatures_updated.csv'  # Provide the path to your CSV file

    data = []
    with open(csv_file, 'r') as file:
        csv_reader = csv.DictReader(file)
        for row in csv_reader:
            selected_row = {col: row[col] for col in selected_columns}
            data.append(selected_row)

    return render(request, 'results_vacsolml.html', {'data': data})
        
def download_csv(request):
    if request.method == 'POST':
        selected_files = request.POST.getlist('csv_files')
        if selected_files:
            # Set the response content type as a zip file
            response = HttpResponse(content_type='application/zip')
            response['Content-Disposition'] = 'attachment; filename="analysis_results.zip"'
            # Create a zip file to store the selected CSV files
            with zipfile.ZipFile(response, 'w') as zip_file:
                # Path to the directory containing the CSV files
                script_dir = sys.path[0]
                directory = os.path.join(script_dir, 'eskapeml/Analysis_Results')
                for file_name in selected_files:
                    # Path to the CSV file
                    file_path = os.path.join(directory, file_name)
                    # Add the CSV file to the zip file
                    zip_file.write(file_path, os.path.basename(file_path))
            return response
    else:
        return render(request, 'results_vacsolml.html')
    
def delete_files(request):
    script_dir = sys.path[0]
    folder_path = os.path.join(script_dir, 'eskapeml/Temp_Results')

    # Delete all files in the folder
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        if os.path.isfile(file_path):
            os.remove(file_path)

    return HttpResponseRedirect('/')    