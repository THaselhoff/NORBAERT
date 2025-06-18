#Python 3.11.7
import threading
import os
import sys
import time
import tkinter as tk #8.6
from tkinter import filedialog #8.6
from tkinter import messagebox #8.6
from tkinter import ttk #8.6
import numpy #1.26.4
import scipy.io.wavfile as wavfile #1.11.4
import matplotlib
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt #3.8.0
matplotlib.use('Agg')
import pandas as pd #2.1.4
from pathlib import Path #1.0.1
import networkx as nx #3.1


#######################################################################

#Functions for Text on input fields while hovering with the cursor
class ToolTip:
    def __init__(self, widget, text):
        self.widget = widget
        self.text = text
        self.tooltip_window = None
        widget.bind("<Enter>", self.show_tooltip)
        widget.bind("<Leave>", self.hide_tooltip)

    def show_tooltip(self, event):
        if self.tooltip_window or not self.text:
            return
        x, y, _, _ = self.widget.bbox("insert")
        x += self.widget.winfo_rootx() + 25
        y += self.widget.winfo_rooty() + 25
        self.tooltip_window = tw = tk.Toplevel(self.widget)
        tw.wm_overrideredirect(True)
        tw.wm_geometry(f"+{x}+{y}")
        label = tk.Label(tw, text=self.text, background="white", relief="solid", borderwidth=1)
        label.pack()

    def hide_tooltip(self, event):
        if self.tooltip_window:
            self.tooltip_window.destroy()
        self.tooltip_window = None

##########################################################################
##########################################################################
##########################################################################
stop_flag = threading.Event()

def start_calculation_thread():
    thread = threading.Thread(target=start_calculation, args=(stop_flag,), daemon=True)
    thread.start()

def close_window():
    stop_flag.set()
    root.destroy()
    os._exit(0)
    
def stop_calculation():
    stop_flag.set()

def browse_directory(var):
    directory = filedialog.askdirectory()
    var.set(directory)
    
def create_AdjMatrix(CorMatrix,threshold):
    cor_matrix=numpy.array(CorMatrix)
    numpy.fill_diagonal(cor_matrix, 0)
    cor_matrix=numpy.where(cor_matrix>=threshold,cor_matrix,0)
    cor_matrix=numpy.where(cor_matrix<threshold,cor_matrix,1)
    adj_matrix=cor_matrix
    adj_matrix=pd.DataFrame(cor_matrix,index=CorMatrix.index,columns=CorMatrix.index)
    return(adj_matrix)

#calculate dB correction from dBFullScale.
def dB_correction(bit_d,dBFS):
    dB_correction=(2**bit_d/2)*10**(-dBFS/20)
    return dB_correction

#calculate root mean square
def RMS(signal):
    rms=numpy.sqrt(numpy.mean(numpy.square(signal)))
    return rms

#Get actual energy level from FFT amplitudes
def RMS_freq(rfft_result, signal):
    rms_freq=numpy.sqrt(numpy.mean(numpy.abs(rfft_result) ** 2/len(signal)))
    return rms_freq

#calculate actual dBZ from RMS, considering the dB_correction
def Lp(messwert, dB_correction):
    Lp=20*numpy.log10(messwert/dB_correction)
    return Lp

def show_warning(message):
    # Create a Tkinter window and hide it (only used to show the dialog)
    root = tk.Tk()
    root.withdraw()
    # Show the warning message
    messagebox.showwarning("Input Error", message)
    # Destroy the root window after the message is closed
    root.destroy()

def convert_seconds(seconds):
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return int(hours), int(minutes), int(seconds)

def start_calculation(stop_flag):
    stop_flag.clear() # Reset the stop flag
     
    input_directory = input_dir_var.get()
    sampling_rate = hz_var.get()
    snippet_length = snippet_length_var.get()
    threshold = threshold_var.get()
    output_directory = output_dir_var.get()
    upper_limit=upper_limit_var.get()

    selected_options = {
        "on_all": on_all_rec.get(),
        "FCM": fcm_var.get(),
        "ADJMatrix": adj_matrix_var.get(),
        "FCM_pdf": fcm_pdf_var.get(),
        "AM_pdf": am_pdf_var.get(),
        "Link Density": link_density_var.get(),

    }
    
    #Warnings
    if len(input_directory)==0:
        show_warning("No input directory is given!")
        return
    
    if len(output_directory)==0:
        show_warning("No output directory is given!")
        return
    
    flist = [os.path.join(dirpath,filename) for dirpath, _, filenames in os.walk(input_directory) for filename in filenames if filename.endswith('.wav')]

    
    if len(flist)==1 and selected_options.get("on_all") is True:
        show_warning("Only one file in directory for between recordings mode!")
        return
    
    if not any([sampling_rate]):
        show_warning("Some input options are empty!")
        return
    
    if not all([upper_limit, snippet_length, threshold]):
        show_warning("Some calculation options are empty!")
        return
    
    if not any(selected_options.values()):
        show_warning("You did not choose an output option!")
        return
    
    if upper_limit > sampling_rate/2:
        show_warning("Your Frequency Limit is higher than the Nyquist Limit! (Sampling Frequency/2)")
        return
    
    if len(flist)==0:
        show_warning("No WAVE files found!")
        return
    

    
    #define constant dbFS for all recordings
    dbfs=129.3 #arbritatrary number, kept for potential future functions which might require correct calibrations
    bit_depth=16 ##arbritatrary number, kept for potential future functions which might require correct calibrations
    #define R2 threshold for FCM to adjacency matrix 
    thresh=threshold**2
  
    #Create Folder paths for output
    if selected_options.get("FCM") is True:
        new_folder_path = Path(output_directory+"/FCMs")
        # Create the directory
        try:
            new_folder_path.mkdir(parents=True, exist_ok=True)  # Creates the directory and any necessary intermediate directories
            print(f"Directory '{new_folder_path}' created successfully.")
        except OSError as e:
            print(f"Error: {e}")   
   
    if selected_options.get("ADJMatrix") is True:
        new_folder_path = Path(output_directory+"/Adj_Matrices")
        # Create the directory
        try:
            new_folder_path.mkdir(parents=True, exist_ok=True)  # Creates the directory and any necessary intermediate directories
            print(f"Directory '{new_folder_path}' created successfully.")
        except OSError as e:
            print(f"Error: {e}") 
  
    
    if selected_options.get("FCM_pdf") is True:
        new_folder_path = Path(output_directory+"/FCM_pdf")
        # Create the directory
        try:
            new_folder_path.mkdir(parents=True, exist_ok=True)  # Creates the directory and any necessary intermediate directories
            print(f"Directory '{new_folder_path}' created successfully.")
        except OSError as e:
            print(f"Error: {e}")


    if selected_options.get("AM_pdf") is True:
        new_folder_path = Path(output_directory+"/Adj_Matrices_pdf")
        # Create the directory
        try:
            new_folder_path.mkdir(parents=True, exist_ok=True)  # Creates the directory and any necessary intermediate directories
            print(f"Directory '{new_folder_path}' created successfully.")
        except OSError as e:
            print(f"Error: {e}")

    #create dataframe to store network measures into
    df_measures=pd.DataFrame(columns=["Link Density"])
    #Progression bar parameters
    prog_bar_counter=1
    total_files = len(flist)
    
    #######################################################      
    #Procedure for calculations between all recordings  
    #######################################################
    if selected_options.get("on_all") is True:
        dB_df_all=pd.DataFrame()
        folder_name=os.path.basename(input_directory) #needed to create final file name
        for file in flist:
            if stop_flag.is_set()==False:                
                try:
                    start_time = time.time()
                    s_rate, sig=wavfile.read(file)
                    try:
                        sig=sig[:,0] #check if mono or stereo file
                    except:
                        sig=sig
                    sig=sig-numpy.mean(sig) #correct for potential dc offsets
                    recording_length=len(sig) 
                    #get frequency labels for each value from rfft
                    freq= numpy.fft.rfftfreq(recording_length, d=1.0/sampling_rate) 
                    #find place where frequency range exceeds 22050 Hz (to stay in line with previous research with following binning)
                    max_freq=numpy.where(freq>=22050)[0][0] 
                    #limit freq labels to 22050 Hz
                    freq=freq[1:max_freq+1] 
                    #Define labels for frequency bins
                    freq_mids=numpy.array_split(freq, 1024)
                    freq_max=[numpy.max(chunks) for chunks in freq_mids]
                    freq_mids=[numpy.median(chunks) for chunks in freq_mids]
                    freq_mids=numpy.array(freq_mids).astype(int)
                    #calculate rfft
                    rfft_result = numpy.fft.rfft(sig)
                    #limit rfft to 22050 kHz
                    rfft_result=rfft_result[1:max_freq] 
                    rfft_bins = numpy.array_split(rfft_result, 1024) #split rfft results into bins 
                    k=dB_correction(bit_depth, dbfs) #calculate dB_correction for given dBFS
                    dB_means = [Lp(RMS_freq(chunk, sig), k) for chunk in rfft_bins] #calculate dB for each bin
                    #store the value at the correct place in dataset place
                    dB_df_all[file.split('\\')[-1]]=dB_means                
                    #update prog bar after each progressed file
                    elapsed_time = time.time() - start_time   
                    
                    execution_time = elapsed_time*(len(flist)-prog_bar_counter) # Calculate the execution time     
                    hours, minutes, seconds = convert_seconds(execution_time)
                    
                    progress_var.set(prog_bar_counter / len(flist) * 100)  # Update progress bar
                    root.update_idletasks()  # Update the GUI
                    counter_label.config(
                        text=f"Files Processed: {prog_bar_counter} / {total_files}\nEstimated Time Remaining: {hours:02d} h {minutes:02d} m {seconds:02d} s",
                        anchor="w", 
                        justify="left" 
                    )

                    prog_bar_counter=prog_bar_counter+1
                   
                    if stop_flag.is_set()==True:
                        break
                    
                except:
                    pass 
         
        if stop_flag.is_set()==False:          
            cor=dB_df_all.T.corr(method="pearson")
            cor.index=freq_mids
            kHz13k=numpy.where(numpy.array(freq_max)>=upper_limit)[0][0]
            AdjMatrix=create_AdjMatrix(cor.iloc[:kHz13k,:kHz13k]**2,thresh)
            if selected_options.get("FCM") is True:
                out=cor.iloc[:kHz13k,:kHz13k]**2
                out.to_csv(output_directory+"/FCMs/"+folder_name+".csv")
            if selected_options.get("ADJMatrix") is True:
                AdjMatrix.astype("int").to_csv(output_directory+"/Adj_Matrices/"+folder_name+".csv")   
                 
            if selected_options.get("FCM_pdf") is True:
                CorMatrix=cor.iloc[:kHz13k,:kHz13k]**2
                ticks=numpy.linspace(0,len(CorMatrix)-1,5).astype("int")
                
                cor_fig=plt.figure()
                ax = cor_fig.add_subplot(111)
                im=ax.imshow(CorMatrix,origin="lower",vmin=0,vmax=1)
                ax.set_title("Frequency Correlation Matrix",fontsize=10, fontweight="bold")
                ax.set_xlabel("Frequency in Hz")
                ax.set_ylabel("Frequency in Hz")
                ax.set_yticks(ticks,labels=CorMatrix.index[ticks])
                ax.set_xticks(ticks,labels=CorMatrix.index[ticks],rotation=45)
                cbar = cor_fig.colorbar(im, ax=ax, label='$R^2$')
                plt.tight_layout()
                cor_fig.savefig(output_directory+"/FCM_pdf/"+folder_name+".pdf", dpi=300)
                plt.close()
                cor_fig.clear()
                del cor_fig
                
            if selected_options.get("AM_pdf") is True:             
                CorMatrix=cor.iloc[:kHz13k,:kHz13k]**2
                cor_matrix=numpy.array(cor.iloc[:kHz13k,:kHz13k])
                numpy.fill_diagonal(cor_matrix, 0)
                cor_matrix=numpy.where(cor_matrix>=threshold,cor_matrix,0)
                cor_matrix=numpy.where(cor_matrix<threshold,cor_matrix,1)
                adj_matrix=cor_matrix
                adj_matrix=pd.DataFrame(cor_matrix,index=CorMatrix.index,columns=CorMatrix.index)
                ticks=numpy.linspace(0,len(cor_matrix)-1,5).astype("int")
                
                adj_fig=plt.figure()
                ax = adj_fig.add_subplot(111)
                ax.imshow(adj_matrix,origin="lower",cmap="binary")
                ax.set_title("Adjacency Matrix (r="+str(threshold)+")",fontsize=10, fontweight="bold")
                ax.set_xlabel("Frequency in Hz")
                ax.set_ylabel("Frequency in Hz")
                ax.set_yticks(ticks,labels=CorMatrix.index[ticks])
                ax.set_xticks(ticks,labels=CorMatrix.index[ticks],rotation=45)
                plt.tight_layout()
                adj_fig.savefig(output_directory+"/Adj_Matrices_pdf/"+folder_name+".pdf", dpi=300)
                plt.close()
                adj_fig.clear()
                del adj_fig
    
            AdjMatrix.index=AdjMatrix.index.astype(str)
            adjmatr=numpy.array(AdjMatrix)
            net=nx.from_numpy_array(adjmatr)
            if selected_options.get("Link Density") is True:
                ld=nx.density(net)
            else:
                ld=numpy.nan

            if selected_options.get("Link Density") is True:
                df_measures.loc[folder_name]=[ld]
                df_measures.to_csv(output_directory+"/Network_Measures_"+folder_name+".csv")
    
            show_warning("Calculation between "+str(prog_bar_counter-1)+" files successfull!")
    
    #######################################################      
    #Procedure for calculations on seperate recordings  
    #######################################################
    elif selected_options.get("on_all") is False: 
    
        for file in flist:
            if stop_flag.is_set()==False:
                start_time = time.time()
                try:
                    s_rate, sig=wavfile.read(file)
                    try:
                        sig=sig[:,0] #check if mono or stereo file
                    except:
                        sig=sig
                    sig=sig-numpy.mean(sig)
                    file_name=os.path.basename(file)
                    recording_length=len(sig)
                    #define maximum recording length that can be divided into 10s clips
                    trunc_recording_len=numpy.floor(recording_length/sampling_rate)*sampling_rate
                    no_10s_clips=(trunc_recording_len/sampling_rate/snippet_length).astype("int")
                    #calculate frequency labels, sampling_rate*snippet length because of X sec clips
                    freq= numpy.fft.rfftfreq(sampling_rate*snippet_length, d=1.0/sampling_rate) 
                    max_freq=numpy.where(freq>=22050)[0][0] 
                    #limit freq labels to 22050 Hz
                    freq=freq[1:max_freq+1] 
                    #Define labels for frequency bins
                    freq_mids=numpy.array_split(freq, 1024)
                    freq_max=[numpy.max(chunks) for chunks in freq_mids]
                    freq_mids=[numpy.median(chunks) for chunks in freq_mids]
                    freq_mids=numpy.array(freq_mids).astype(int)
                except:
                    pass 
                dB_df=pd.DataFrame(index=freq_mids,columns=numpy.arange(no_10s_clips))
                sig_cut=sig[:trunc_recording_len.astype("int")]
                seg=sampling_rate*snippet_length
                for i in range(no_10s_clips):
                    sig_seg=sig_cut[i*seg:i*seg+seg]
                
                    rfft_result = numpy.fft.rfft(sig_seg) #calculate rfft
                    rfft_result=rfft_result[1:max_freq] #limit rfft to 22050
                    rfft_bins = numpy.array_split(rfft_result, 1024) #split rfft results into bins 
                    k=dB_correction(bit_depth, dbfs) #calculate dB_correction for given dBFS
                    dB_means = [Lp(RMS_freq(chunk, sig), k) for chunk in rfft_bins]
                    dB_df[i]=dB_means #store the value at the correct place in dataset place
                else:
                    pass
                
                cor=dB_df.T.corr(method="pearson")
                kHz13k=numpy.where(numpy.array(freq_max)>=upper_limit)[0][0]
                AdjMatrix=create_AdjMatrix(cor.iloc[:kHz13k,:kHz13k]**2,thresh)
                if selected_options.get("FCM") is True:
                    out=cor.iloc[:kHz13k,:kHz13k]**2
                    out.to_csv(output_directory+"/FCMs/"+file_name+".csv")
                if selected_options.get("ADJMatrix") is True:
                    AdjMatrix.to_csv(output_directory+"/Adj_Matrices/"+file_name+".csv")   
                     
                if selected_options.get("FCM_pdf") is True:
                    CorMatrix=cor.iloc[:kHz13k,:kHz13k]**2
                    ticks=numpy.linspace(0,len(CorMatrix)-1,5).astype("int")
                    
                    cor_fig=plt.figure()
                    ax = cor_fig.add_subplot(111)
                    im=ax.imshow(CorMatrix,origin="lower",vmin=0,vmax=1)
                    ax.set_title("Frequency Correlation Matrix",fontsize=10, fontweight="bold")
                    ax.set_xlabel("Frequency in Hz")
                    ax.set_ylabel("Frequency in Hz")
                    ax.set_yticks(ticks,labels=CorMatrix.index[ticks])
                    ax.set_xticks(ticks,labels=CorMatrix.index[ticks],rotation=45)
                    cbar = cor_fig.colorbar(im, ax=ax, label='$R^2$')
                    plt.tight_layout()
                    cor_fig.savefig(output_directory+"/FCM_pdf/"+file_name+".pdf", dpi=300)
                    plt.close()
                    cor_fig.clear()
                    del cor_fig
                    
                if selected_options.get("AM_pdf") is True:             
                    CorMatrix=cor.iloc[:kHz13k,:kHz13k]**2
                    cor_matrix=numpy.array(cor.iloc[:kHz13k,:kHz13k])
                    numpy.fill_diagonal(cor_matrix, 0)
                    cor_matrix=numpy.where(cor_matrix>=threshold,cor_matrix,0)
                    cor_matrix=numpy.where(cor_matrix<threshold,cor_matrix,1)
                    adj_matrix=cor_matrix
                    adj_matrix=pd.DataFrame(cor_matrix,index=CorMatrix.index,columns=CorMatrix.index)
                    ticks=numpy.linspace(0,len(cor_matrix)-1,5).astype("int")
                    
                    adj_fig=plt.figure()
                    ax = adj_fig.add_subplot(111)
                    ax.imshow(adj_matrix,origin="lower",cmap="binary")
                    ax.set_title("Adjacency Matrix (r="+str(threshold)+")",fontsize=10, fontweight="bold")
                    ax.set_xlabel("Frequency in Hz")
                    ax.set_ylabel("Frequency in Hz")
                    ax.set_yticks(ticks,labels=CorMatrix.index[ticks])
                    ax.set_xticks(ticks,labels=CorMatrix.index[ticks],rotation=45)
                    plt.tight_layout()
                    adj_fig.savefig(output_directory+"/Adj_Matrices_pdf/"+file_name+".pdf", dpi=300)
                    plt.close()
                    adj_fig.clear()
                    del adj_fig
        
                AdjMatrix.index=AdjMatrix.index.astype(str)
                adjmatr=numpy.array(AdjMatrix)
                net=nx.from_numpy_array(adjmatr)
                if selected_options.get("Link Density") is True:
                    ld=nx.density(net)
                else:
                    ld=numpy.nan
                df_measures.loc[file_name]=[ld]
                
                #update prog bar
                elapsed_time = time.time() - start_time   
                
                execution_time = elapsed_time*(len(flist)-prog_bar_counter) # Calculate the execution time     
                hours, minutes, seconds = convert_seconds(execution_time)
                
                progress_var.set(prog_bar_counter / len(flist) * 100)  # Update progress bar
                root.update_idletasks()  # Update the GUI
                counter_label.config(
                    text=f"Files Processed: {prog_bar_counter} / {total_files}\nEstimated Time Remaining: {hours:02d} h {minutes:02d} m {seconds:02d} s",
                    anchor="w", 
                    justify="left"  
                )


                prog_bar_counter=prog_bar_counter+1
            else:
                pass
            if selected_options.get("Link Density") is True:
                df_measures.to_csv(output_directory+"/Network_Measures.csv")
            if stop_flag.is_set()==True:
                break
        show_warning("Calculation of "+str(prog_bar_counter-1)+" files successfull!")

#######################
# UI
#######################   
# Create the main window
root = tk.Tk()
root.title("NORBAERT")

style = ttk.Style()

# Create variables to hold input data
input_dir_var = tk.StringVar()
hz_var = tk.IntVar()
recording_length_var = tk.IntVar()
upper_limit_var=tk.DoubleVar()
snippet_length_var = tk.IntVar()
threshold_var = tk.DoubleVar()
output_dir_var = tk.StringVar()

on_all_rec= tk.BooleanVar()
fcm_var = tk.BooleanVar()
adj_matrix_var = tk.BooleanVar()
fcm_pdf_var = tk.BooleanVar()
am_pdf_var = tk.BooleanVar()
link_density_var = tk.BooleanVar()


#######################
# Directories Section
#######################
dir_frame_label = ttk.Label(root, text="Directories", font=('Helvetica', 8, 'bold'))
dir_frame = ttk.LabelFrame(root, labelwidget=dir_frame_label)
dir_frame.pack(pady=5, padx=10, fill='x')

tk.Label(dir_frame, text="Input Directory:").grid(row=0, column=0, sticky='w', padx=5)
input_dir_entry = ttk.Entry(dir_frame, textvariable=input_dir_var, width=30)
input_dir_entry.grid(row=0, column=1)
ttk.Button(dir_frame, text="Browse", command=lambda: browse_directory(input_dir_var)).grid(row=0, column=2)
ToolTip(input_dir_entry, "Directory containing the WAVE files to process.")

tk.Label(dir_frame, text="Output Directory:").grid(row=1, column=0, sticky='w',padx=5)
output_dir_entry = ttk.Entry(dir_frame, textvariable=output_dir_var, width=30)
output_dir_entry.grid(row=1, column=1)
ttk.Button(dir_frame, text="Browse", command=lambda: browse_directory(output_dir_var)).grid(row=1, column=2)
ToolTip(output_dir_entry, "Directory where output files will be saved.")

all_rec = tk.Checkbutton(dir_frame, text="Between all recordings", bg="lightgrey", variable=on_all_rec)
all_rec.grid(row=2, column=0, columnspan=3, sticky='w', padx=5, pady=(10,10))
ToolTip(all_rec, "Calculate one FCM and Adjacency Matrix between all recordings in specified input directory.\nIf checked, this will not produce outputs for single recordings.")

#######################
# Input Options
#######################
input_frame_label = ttk.Label(root, text="Input Options", font=('Helvetica', 8, 'bold'))
input_frame = ttk.LabelFrame(root, labelwidget=input_frame_label)
input_frame.pack(pady=5, padx=10, fill='x')

tk.Label(input_frame, text="All recordings need to have the same sampling frequency.", fg='dimgray',font=('Helvetica',8, 'italic'), justify='left').grid(row=0, column=0, columnspan=3, sticky='w', padx=5)
#
hz_var.set(44100)
tk.Label(input_frame, text="Sampling Frequency:").grid(row=1, column=0, sticky='w', padx=5)
#
hz_entry = ttk.Entry(input_frame, textvariable=hz_var, width=6, justify=tk.CENTER)
tk.Label(input_frame, text="Hz").grid(row=1, column=2, sticky='w', pady=(0,5))
hz_entry.grid(row=1, column=1)
ToolTip(hz_entry, "Sampling rate of the audio files in Hz (min. 44,100).")

#######################
# Calculation Parameters
#######################
calc_frame_label = ttk.Label(root, text="Calculation Parameters", font=('Helvetica', 8, 'bold'))
calc_frame = ttk.LabelFrame(root, labelwidget=calc_frame_label)
calc_frame.pack(pady=5, padx=10, fill='x')

upper_limit_var.set(13000)  # Preset the upper freq limit to 13 kHz
tk.Label(calc_frame, text="Frequency Limit:").grid(row=0, column=0, sticky='w', padx=5)
uper_limit_entry = ttk.Entry(calc_frame, textvariable=upper_limit_var, width=6, justify=tk.CENTER)
tk.Label(calc_frame, text="Hz").grid(row=0, column=2, sticky='w')
uper_limit_entry.grid(row=0, column=1)
ToolTip(uper_limit_entry, "Upper Frequency limit for the FCM.")

snippet_length_var.set(10)  # Preset the snippet length to 10
tk.Label(calc_frame, text="Snippet Length:").grid(row=1, column=0, sticky='w', padx=5)
snippet_length_entry = ttk.Entry(calc_frame, textvariable=snippet_length_var, width=6, justify=tk.CENTER)
tk.Label(calc_frame, text="s").grid(row=1, column=2, sticky='w')
snippet_length_entry.grid(row=1, column=1)
ToolTip(snippet_length_entry, "Length of each snippet in seconds.")

threshold_var.set(0.8)  # Preset the snippet length to 10
tk.Label(calc_frame, text="Threshold:").grid(row=2, column=0, sticky='w', padx=5)
threshold_entry = ttk.Entry(calc_frame, textvariable=threshold_var, width=6, justify=tk.CENTER)
tk.Label(calc_frame, text="r").grid(row=2, column=2, sticky='w')
threshold_entry.grid(row=2, column=1, pady=(0,5))
ToolTip(threshold_entry, "Pearson correlation threshold for the Frequency Correlation Matrix to obtain the Adjacency Matrix.")

#######################
# Output Options
#######################
output_frame_label = ttk.Label(root, text="Output Options", font=('Helvetica', 8, 'bold'))
output_frame = ttk.LabelFrame(root, labelwidget=output_frame_label)
output_frame.pack(pady=5, padx=10, fill='x')

#
tk.Label(output_frame, text="Tables (.csv)", fg='dimgray',font=('Helvetica',8, 'italic')).grid(row=1, column=0, sticky='w', padx=5)
#
fcm_check = tk.Checkbutton(output_frame, text="FCM", variable=fcm_var)
fcm_check.grid(row=2, column=0, sticky='w', padx=5)
ToolTip(fcm_check, "Include Frequency Correlation Matrix in the output as .csv.")

adj_matrix_check = tk.Checkbutton(output_frame, text="Adjacency Matrix", variable=adj_matrix_var)
adj_matrix_check.grid(row=2, column=1, sticky='w', padx=5)
ToolTip(adj_matrix_check, "Include Adjacency Matrix in the output as .csv.")

#
tk.Label(output_frame, text="Network measures (.csv)", fg='dimgray',font=('Helvetica',8, 'italic')).grid(row=3, column=0, sticky='w', padx=5)
#
link_density_check = tk.Checkbutton(output_frame, text="Link Density", variable=link_density_var)
link_density_check.grid(row=4, column=0, sticky='w', padx=5)
ToolTip(link_density_check, "Calculate Link Density.")


#
tk.Label(output_frame, text="Figures (.pdf)", fg='dimgray',font=('Helvetica',8, 'italic')).grid(row=5, column=0, sticky='w', padx=5)
#
fcm_check = tk.Checkbutton(output_frame, text="FCM", variable=fcm_pdf_var)
fcm_check.grid(row=6, column=0, sticky='w', padx=5)
ToolTip(fcm_check, "Include Frequency Correlation Matrix in the output. Note that this may require significant hard drive space.")

adj_matrix_check = tk.Checkbutton(output_frame, text="Adjacency Matrix", variable=am_pdf_var)
adj_matrix_check.grid(row=6, column=1, sticky='w', padx=5)
ToolTip(adj_matrix_check, "Include Adjacency Matrix in the output. Note that this may require significant hard drive space.")

#######################
# Control Buttons
#######################

control_frame = ttk.Frame(root)
control_frame.pack(pady=10, padx=10, fill='x')

ttk.Button(control_frame, text="Start Calculation", command=start_calculation_thread).grid(row=0, column=2, pady=(15, 0), padx=(5,5))
ttk.Button(control_frame, text="Stop Calculation", command=stop_calculation).grid(row=1, column=2, pady=(0,15), padx=(5,5))
ttk.Button(control_frame, text="Close", command=close_window).grid(row=2, column=0, columnspan=3, pady=(0,15))

root.update_idletasks()

progress_var = tk.DoubleVar()
progress_bar = ttk.Progressbar(control_frame, length=260, variable=progress_var, maximum=100)
progress_bar.grid(row=0, column=0, columnspan=2, padx=(5,0), pady=(15, 0), sticky="w")

counter_label = tk.Label(control_frame, text="Files Processed: 0 / 0")
counter_label.grid(row=1, column=0, columnspan=2, padx=(5,0), pady=(0,15), sticky="w")

tk.Label(control_frame, text="v. 0.1.0", fg='dimgray',font=('Helvetica',8, 'italic')).grid(row=2, column=2, sticky='e', padx=5, pady=(20,0))



##################################
root.update_idletasks()

w = root.winfo_width() # width for the Tk root
h = root.winfo_height() # height for the Tk root
# get screen width and height
ws = root.winfo_screenwidth() # width of the screen
hs = root.winfo_screenheight() # height of the screen
# calculate x and y coordinates for the Tk root window
x = (ws/2) - (w/2)
y = (hs/2) - (h/2)

root.geometry('%dx%d+%d+%d' % (w, h, x, y))
# Run the Tkinter event loop
root.mainloop()


