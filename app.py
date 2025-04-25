#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python
# coding: utf-8

import os
import io
import re
import textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import xlsxwriter
from flask import Flask, request, render_template, redirect, url_for, flash, send_from_directory
from werkzeug.utils import secure_filename
from scipy.stats import ttest_ind, f_oneway

app = Flask(__name__)
app.config['SECRET_KEY'] = 'your-secret-key'
app.config['UPLOAD_FOLDER'] = 'uploads'
app.config['STATIC_FOLDER'] = 'static'
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
os.makedirs(app.config['STATIC_FOLDER'], exist_ok=True)

# --- Helper function to sanitize Excel sheet names ---
def sanitize_sheetname(name):
    """
    Replace invalid Excel sheet name characters with an underscore and limit to 31 characters.
    Invalid characters include: : \ / ? * [ ]
    """
    sanitized = re.sub(r'[:\\/?*\[\]]', '_', name)
    return sanitized[:31]

# --- Helper function to annotate significance on a bar graph ---
def annotate_significance(ax, x1, x2, y, h, significance_text):
    """Draws a horizontal line between x1 and x2 at height y+h and adds significance text above it."""
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c='k')
    ax.text((x1+x2)*.5, y+h, significance_text, ha='center', va='bottom', color='k', fontsize=10)

# --- Helper function to compute ranking and pairwise p-values for a given metric ---
def compute_ranking_and_pvalues(metric_values):
    """
    metric_values: dict mapping sample (string) to a list of replicate values.
    Returns:
      ranking: dict mapping rank number (starting at 1) to list of sample names that are not significantly different.
      pvals: dict with keys (sample1, sample2) and p-value.
    """
    # Compute mean for each sample
    means = {sample: np.nanmean(vals) for sample, vals in metric_values.items()}
    # Sort samples by mean (descending)
    sorted_samples = sorted(means.items(), key=lambda x: x[1], reverse=True)
    groups = []
    for sample, _ in sorted_samples:
        placed = False
        for group in groups:
            rep = group[0]
            stat, p = ttest_ind(metric_values[sample], metric_values[rep], equal_var=False, nan_policy='omit')
            if p >= 0.05:
                group.append(sample)
                placed = True
                break
        if not placed:
            groups.append([sample])
    ranking = {}
    for i, group in enumerate(groups, start=1):
        ranking[i] = group
    # Compute all pairwise p-values
    samples = list(metric_values.keys())
    pvals = {}
    for i in range(len(samples)):
        for j in range(i+1, len(samples)):
            s1 = samples[i]
            s2 = samples[j]
            stat, p = ttest_ind(metric_values[s1], metric_values[s2], equal_var=False, nan_policy='omit')
            pvals[(s1, s2)] = float(f"{p:.8f}")
    return ranking, pvals

def process_file(file, effective_electrode, params, electrode_actual, day_value):
    # Save uploaded file securely
    filename = secure_filename(file.filename)
    file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
    file.save(file_path)
    
    # Read file (CSV or XLSX)
    try:
        if filename.lower().endswith('.csv'):
            df_raw = pd.read_csv(file_path, header=None)
        else:
            df_raw = pd.read_excel(file_path, header=None)
    except Exception as e:
        raise Exception(f"Error reading file {filename}: {str(e)}")
    
    # --- Process Raw Data (assume data starts at row 4) ---
    data = df_raw.iloc[3:].copy()
    data.columns = ['raw']
    def split_data(cell):
        try:
            parts = str(cell).split(';')
            return pd.Series([float(parts[0]), float(parts[1].replace('"',''))])
        except:
            return pd.Series([np.nan, np.nan])
    data[['Time/s', 'Current/uA']] = data['raw'].apply(split_data)
    data = data.dropna(subset=['Time/s']).reset_index(drop=True)
    df_data = data[['Time/s', 'Current/uA']]
    
    # --- Generate Sampling Timestamps ---
    sampling_times = [params['time_glu1'] - 50]
    total_drops = params['num_drops_glu1'] + params['num_drops_glu2']
    for drop in range(1, total_drops+1):
        t = params['time_glu1'] + params['drop_interval']/2 + (effective_electrode-1)*params['time_delay'] + params['drop_interval']*(drop-1)
        sampling_times.append(t)
    df_sampling = pd.DataFrame({'Time Sampled/s': sampling_times})
    
    # --- Get Current Values at Sampling Times ---
    sampled_current = []
    for t in sampling_times:
        idx = (np.abs(df_data['Time/s'] - t)).idxmin()
        sampled_current.append(df_data.loc[idx, 'Current/uA'])
    df_sampling['Current/uA'] = sampled_current
    
    # --- Compute Glucose Addition & Final Concentration ---
    vol1_list, vol2_list, amount_list, tot_vol_list, final_conc_list = [], [], [], [], []
    for i in range(len(df_sampling)):
        if i == 0:
            vol1 = 0
            vol2 = 0
        else:
            drop_index = i
            vol1 = params['vol_glu1'] * drop_index if drop_index <= params['num_drops_glu1'] else params['vol_glu1'] * params['num_drops_glu1']
            vol2 = params['vol_glu2'] * (drop_index - params['num_drops_glu1']) if drop_index > params['num_drops_glu1'] else 0
        vol1_list.append(vol1)
        vol2_list.append(vol2)
        amt = vol1 * params['conc_glu1'] + vol2 * params['conc_glu2']
        amount_list.append(amt)
        tot_vol = params['vol_sweat'] + vol1 + vol2
        tot_vol_list.append(tot_vol)
        final_conc = amt / tot_vol if tot_vol != 0 else 0
        final_conc_list.append(final_conc)
    col_vol1 = f"Volume of {params['conc_glu1']}mM drops/uL"
    col_vol2 = f"Volume of {params['conc_glu2']}mM drops/uL"
    df_sampling[col_vol1] = vol1_list
    df_sampling[col_vol2] = vol2_list
    df_sampling['Amount of Glucose Added/nmol'] = amount_list
    df_sampling['Total Volume/uL'] = tot_vol_list
    df_sampling['Final Conc. Of Glucose in the Solution on the Electrode/mM'] = final_conc_list
    df_sampling.rename(columns={
        'Current/uA': 'Current at Time Sampled/uA',
        col_vol1: f"Volume of {params['conc_glu1']}mM Glucose Solution in the Droplet on the Electrode/uL",
        col_vol2: f"Volume of {params['conc_glu2']}mM Glucose Solution in the Droplet on the Electrode/uL",
        'Amount of Glucose Added/nmol': "Amount of Glucose in the Droplet on the Electrode/nmol",
        'Total Volume/uL': "Total Volume of the Droplet on the Electrode/uL",
        'Final Conc. Of Glucose in the Solution on the Electrode/mM': "Concentration of Glucose in the Droplet on the Electrode/mM"
    }, inplace=True)
    
    # --- Processed Data: Trendline (linear regression) ---
    x_proc = np.array(final_conc_list)
    y_proc = np.array(df_sampling['Current at Time Sampled/uA'])
    valid = ~np.isnan(y_proc)
    x_proc = x_proc[valid]
    y_proc = y_proc[valid]
    if len(x_proc) > 1:
        coeffs = np.polyfit(x_proc, y_proc, 1)
        slope, intercept = coeffs
        y_fit = np.polyval(coeffs, x_proc)
        ss_res = np.sum((y_proc - y_fit)**2)
        ss_tot = np.sum((y_proc - np.mean(y_proc))**2)
        r2 = 1 - ss_res/ss_tot if ss_tot != 0 else 0
    else:
        slope, intercept, r2 = 0, 0, 0
    sensitivity = abs(slope)
    trendline_text = f"Trendline: y = {slope:.4f}x + {intercept:.4f}, R² = {r2:.4f}"
    sensitivity_text = f"The sensitivity of Electrode {electrode_actual} casted with {params['electrode_cast']} on Day {day_value} is {sensitivity:.4f} mA/M."
    
    # --- Graph Titles ---
    raw_title = f"Graph of Current/uA vs Time/s for Electrode {electrode_actual} casted with {params['electrode_cast']} on Day {day_value}"
    processed_title = f"Graph of Current/uA vs Glucose Concentration/mM\nfor Electrode {electrode_actual} casted with {params['electrode_cast']} on Day {day_value}"
    
    # --- Generate Graphs ---
    # Raw Data Graph
    raw_graph_path = os.path.join(app.config['STATIC_FOLDER'], f"raw_graph_{electrode_actual}_{day_value}.png")
    plt.figure()
    plt.plot(df_data['Time/s'], df_data['Current/uA'], color='blue')
    plt.title(raw_title, wrap=True)
    plt.xlabel("Time/s")
    plt.ylabel("Current/uA")
    if params['x_min'] and params['x_max']:
        plt.xlim(float(params['x_min']), float(params['x_max']))
    if params['y_min'] and params['y_max']:
        plt.ylim(float(params['y_min']), float(params['y_max']))
    plt.grid(True)
    plt.savefig(raw_graph_path, bbox_inches="tight")
    plt.close()
    
    # Raw Data Graph from 500s
    raw_graph_500_path = os.path.join(app.config['STATIC_FOLDER'], f"raw_graph_500_{electrode_actual}_{day_value}.png")
    df_data_500 = df_data[df_data['Time/s'] >= 500]
    plt.figure()
    plt.plot(df_data_500['Time/s'], df_data_500['Current/uA'], color='blue')
    plt.title(raw_title + " (from 500s)", wrap=True)
    plt.xlabel("Time/s")
    plt.ylabel("Current/uA")
    plt.grid(True)
    plt.savefig(raw_graph_500_path, bbox_inches="tight")
    plt.close()
    
    # Processed Data Graph (scatter with trendline)
    processed_graph_path = os.path.join(app.config['STATIC_FOLDER'], f"processed_graph_{electrode_actual}_{day_value}.png")
    plt.figure()
    plt.scatter(df_sampling["Concentration of Glucose in the Droplet on the Electrode/mM"],
                df_sampling["Current at Time Sampled/uA"], color='green')
    if len(x_proc) > 1:
        coeffs2 = np.polyfit(x_proc, y_proc, 1)
        slope2, intercept2 = coeffs2
        x_sorted = np.sort(x_proc)
        y_fit2 = slope2 * x_sorted + intercept2
        plt.plot(x_sorted, y_fit2, 'r--', label=f"y = {slope2:.4f}x + {intercept2:.4f}, R² = {r2:.4f}")
        plt.legend()
    plt.title(processed_title, wrap=True)
    plt.xlabel("Concentration (mM)")
    plt.ylabel("Current (uA)")
    plt.grid(True)
    plt.savefig(processed_graph_path, bbox_inches="tight")
    plt.close()
    
    # --- Generate Excel File for this electrode/day ---
    excel_chart_path = os.path.join(app.config['STATIC_FOLDER'], f"excel_charts_{electrode_actual}_{day_value}.xlsx")
    workbook = xlsxwriter.Workbook(excel_chart_path)
    # Raw Data sheet
    ws_raw = workbook.add_worksheet("Raw Data")
    headers_raw = df_data.columns.tolist()
    ws_raw.write_row(0, 0, headers_raw)
    for i, row in df_data.iterrows():
        ws_raw.write_row(i+1, 0, row.tolist())
    chart1 = workbook.add_chart({'type': 'line'})
    n = len(df_data)
    chart1.add_series({
        'name': 'Raw Data',
        'categories': ['Raw Data', 1, 0, n, 0],
        'values': ['Raw Data', 1, 1, n, 1],
    })
    chart1.set_title({'name': raw_title})
    chart1.set_x_axis({'name': "Time/s"})
    chart1.set_y_axis({'name': "Current/uA"})
    if params['x_min'] and params['x_max']:
        chart1.set_x_axis({'min': float(params['x_min']), 'max': float(params['x_max'])})
    if params['y_min'] and params['y_max']:
        chart1.set_y_axis({'min': float(params['y_min']), 'max': float(params['y_max'])})
    ws_raw.insert_chart('D2', chart1)
    # Processed Data sheet
    ws_proc = workbook.add_worksheet("Processed Data")
    headers_proc = df_sampling.columns.tolist()
    ws_proc.write_row(0, 0, headers_proc)
    for i, row in df_sampling.iterrows():
        ws_proc.write_row(i+1, 0, row.tolist())
    chart2 = workbook.add_chart({'type': 'scatter'})
    n2 = len(df_sampling)
    cat_col = headers_proc.index("Concentration of Glucose in the Droplet on the Electrode/mM")
    val_col = headers_proc.index("Current at Time Sampled/uA")
    chart2.add_series({
        'name': 'Processed Data',
        'categories': ['Processed Data', 1, cat_col, n2, cat_col],
        'values': ['Processed Data', 1, val_col, n2, val_col],
        'trendline': {'type': 'linear', 'display_equation': True, 'display_r_squared': True},
    })
    chart2.set_title({'name': processed_title})
    chart2.set_x_axis({'name': "Concentration (mM)"})
    chart2.set_y_axis({'name': "Current (uA)"})
    ws_proc.insert_chart('H2', chart2)
    workbook.close()
    
    table_html = df_sampling.to_html(classes='data', border=0)
    
    return {
        'electrode': electrode_actual,
        'day': day_value,
        'raw_title': raw_title,
        'processed_title': processed_title,
        'trendline_text': trendline_text,
        'sensitivity_text': sensitivity_text,
        'sensitivity': sensitivity,
        'raw_graph': os.path.basename(raw_graph_path),
        'raw_graph_500': os.path.basename(raw_graph_500_path),
        'processed_graph': os.path.basename(processed_graph_path),
        'excel_chart': os.path.basename(excel_chart_path),
        'table_html': table_html
    }

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        # --- Read Common Parameters ---
        try:
            num_electrodes = int(request.form.get('num_electrodes', 1))
            num_days = int(request.form.get('num_days', 1))
            time_delay = float(request.form.get('time_delay', 4))
            vol_glu1 = float(request.form.get('vol_glu1', 2))
            conc_glu1 = float(request.form.get('conc_glu1', 1))
            time_glu1 = float(request.form.get('time_glu1', 1000))
            vol_glu2 = float(request.form.get('vol_glu2', 2))
            conc_glu2 = float(request.form.get('conc_glu2', 10))
            time_glu2 = float(request.form.get('time_glu2', 2000))
            drop_interval = float(request.form.get('drop_interval', 100))
            num_drops_glu1 = int(request.form.get('num_drops_glu1', 10))
            num_drops_glu2 = int(request.form.get('num_drops_glu2', 10))
            vol_sweat = float(request.form.get('vol_sweat', 50))
            x_min = request.form.get('x_min')
            x_max = request.form.get('x_max')
            y_min = request.form.get('y_min')
            y_max = request.form.get('y_max')
        except Exception as e:
            flash("Error reading common parameters: " + str(e))
            return redirect(request.url)
        params = {
            'time_delay': time_delay,
            'vol_glu1': vol_glu1,
            'conc_glu1': conc_glu1,
            'time_glu1': time_glu1,
            'vol_glu2': vol_glu2,
            'conc_glu2': conc_glu2,
            'time_glu2': time_glu2,
            'drop_interval': drop_interval,
            'num_drops_glu1': num_drops_glu1,
            'num_drops_glu2': num_drops_glu2,
            'vol_sweat': vol_sweat,
            'x_min': x_min,
            'x_max': x_max,
            'y_min': y_min,
            'y_max': y_max
        }
        # --- Read Electrode Info ---
        electrode_infos = []
        for i in range(1, num_electrodes+1):
            try:
                elec_num = int(request.form.get('electrode_number_' + str(i)))
            except:
                elec_num = i
            elec_cast = request.form.get('electrode_cast_' + str(i))
            electrode_infos.append({'actual': elec_num, 'cast': elec_cast})
        # --- Read Day Info ---
        days_list = []
        for j in range(1, num_days+1):
            day_val = request.form.get('day_' + str(j))
            days_list.append(day_val)
        initial_day = days_list[0]
        final_day = days_list[-1]
        
        individual_results = {}
        sensitivity_dict = {}  # electrode -> { day: sensitivity }
        for elec in electrode_infos:
            sensitivity_dict[elec['actual']] = {}
            for day in days_list:
                file_key = f"file_{elec['actual']}_{day}"
                file = request.files.get(file_key)
                if not file:
                    flash(f"File for Electrode {elec['actual']} Day {day} not uploaded.")
                    return redirect(request.url)
                effective = ((elec['actual']-1)%8)+1
                params['electrode_cast'] = elec['cast']
                try:
                    result = process_file(file, effective, params, elec['actual'], day)
                except Exception as e:
                    flash(f"Error processing file for Electrode {elec['actual']} Day {day}: {str(e)}")
                    return redirect(request.url)
                individual_results[(elec['actual'], day)] = result
                sensitivity_dict[elec['actual']][day] = result['sensitivity']
        
        # --- Build Individual Sensitivity Summary Table ---
        summary_rows = []
        for elec in electrode_infos:
            row = {"Electrode": f"Electrode {elec['actual']} casted with {elec['cast']}"}
            sens_values = []
            for day in days_list:
                val = sensitivity_dict[elec['actual']].get(day, np.nan)
                row[f"Sensitivity on Day {day} (mA/M)"] = val
                sens_values.append(val)
            row["Average Sensitivity Across Days (mA/M)"] = np.nanmean(sens_values)
            row[f"Change in Sensitivity (Sensitivity on Day {final_day} - Sensitivity on Day {initial_day}) (mA/M)"] = sensitivity_dict[elec['actual']].get(final_day, np.nan) - sensitivity_dict[elec['actual']].get(initial_day, np.nan)
            summary_rows.append(row)
        summary_df = pd.DataFrame(summary_rows)
        
        # --- Build Aggregated Sensitivity Summary (by Sample) ---
        agg_dict = {}
        sample_sens_all = {}
        for elec in electrode_infos:
            cast = elec['cast']
            if cast not in agg_dict:
                agg_dict[cast] = {day: [] for day in days_list}
                sample_sens_all[cast] = []
            for day in days_list:
                val = sensitivity_dict[elec['actual']].get(day, np.nan)
                agg_dict[cast][day].append(val)
            for idx, row in summary_df.iterrows():
                if elec['cast'] in row["Electrode"]:
                    sample_sens_all[cast].append(row["Average Sensitivity Across Days (mA/M)"])
                    break
        agg_rows = []
        for cast, day_vals in agg_dict.items():
            row = {"Sample Casted on Electrode": cast}
            for day in days_list:
                avg_val = np.nanmean(day_vals[day]) if len(day_vals[day])>0 else np.nan
                row[f"Sensitivity on Day {day} (mA/M)"] = avg_val
            row["Average Sensitivity Across Days (mA/M)"] = np.nanmean([np.nanmean(day_vals[day]) for day in days_list])
            row[f"Change in Sensitivity (Sensitivity on Day {final_day} - Sensitivity on Day {initial_day}) (mA/M)"] = row.get(f"Sensitivity on Day {final_day} (mA/M)", np.nan) - row.get(f"Sensitivity on Day {initial_day} (mA/M)", np.nan)
            agg_rows.append(row)
        aggregated_df = pd.DataFrame(agg_rows)
        
        # --- Build Standard Deviation Table for Aggregated Data ---
        std_rows = []
        for cast, day_vals in agg_dict.items():
            row = {"Sample Casted on Electrode": cast}
            for day in days_list:
                std_val = np.std(day_vals[day], ddof=1) if len(day_vals[day])>1 else 0
                row[f"Standard Deviation of Sensitivity on Day {day} (mA/M)"] = std_val
            row["Standard Deviation in Average Sensitivity Across Days (mA/M)"] = np.std(sample_sens_all[cast], ddof=1) if len(sample_sens_all[cast])>1 else 0
            changes = []
            for elec in electrode_infos:
                if elec['cast'] == cast:
                    diff = sensitivity_dict[elec['actual']].get(final_day, np.nan) - sensitivity_dict[elec['actual']].get(initial_day, np.nan)
                    changes.append(diff)
            row[f"Standard Deviation in Change in Sensitivity (Sensitivity on Day {final_day} - Sensitivity on Day {initial_day}) (mA/M)"] = np.std(changes, ddof=1) if len(changes)>1 else 0
            std_rows.append(row)
        std_df = pd.DataFrame(std_rows)
        
        # --- Compute Rankings and Pairwise P-values for All Sensitivity Metrics ---
        metrics_data = {}
        metrics_data["Average Sensitivity Across Days (mA/M)"] = sample_sens_all
        for day in days_list:
            day_metric = {}
            for cast in agg_dict:
                day_metric[cast] = agg_dict[cast][day]
            metrics_data[f"Sensitivity on Day {day} (mA/M)"] = day_metric
        change_metric = {}
        for cast in agg_dict:
            diffs = []
            for elec in electrode_infos:
                if elec['cast'] == cast:
                    diff = sensitivity_dict[elec['actual']].get(final_day, np.nan) - sensitivity_dict[elec['actual']].get(initial_day, np.nan)
                    diffs.append(diff)
            change_metric[cast] = diffs
        metrics_data[f"Change in Sensitivity (Sensitivity on Day {final_day} - Sensitivity on Day {initial_day}) (mA/M)"] = change_metric
        
        ranking_all = {}
        pairwise_all = {}
        for metric, data_dict in metrics_data.items():
            ranking, pvals = compute_ranking_and_pvalues(data_dict)
            ranking_all[metric] = ranking
            pairwise_all[metric] = pvals
        
        ranking_table_rows = []
        max_groups = max(len(ranks) for ranks in ranking_all.values())
        for metric, ranks in ranking_all.items():
            row = {"Metric": metric}
            for i in range(1, max_groups+1):
                row[f"Rank {i}"] = ", ".join(ranks[i]) if i in ranks else "-"
            ranking_table_rows.append(row)
        ranking_all_df = pd.DataFrame(ranking_table_rows)
        ranking_table_html = ranking_all_df.to_html(classes='data', border=0, index=False)
        
        pvalues_tables_html = {}
        for metric, pvals in pairwise_all.items():
            rows = []
            for (s1, s2), p in pvals.items():
                rows.append({"Sample 1": s1, "Sample 2": s2, "p-value": p})
            df_p = pd.DataFrame(rows)
            pvalues_tables_html[metric] = df_p.to_html(classes='data', border=0, index=False)
        
        # --- Generate Excel Summary for Aggregated Data with extra tabs ---
        combined_excel_path = os.path.join(app.config['STATIC_FOLDER'], "Summary_of_Electrochemical_Sensitivity_of_Samples.xlsx")
        with pd.ExcelWriter(combined_excel_path, engine='xlsxwriter') as writer:
            summary_df.to_excel(writer, sheet_name="Individual Sensitivity Summary", index=False)
            aggregated_df.to_excel(writer, sheet_name="Aggregated Sensitivity Summary", index=False)
            std_df.to_excel(writer, sheet_name="Aggregated Sensitivity StdDev", index=False)
            ranking_all_df.to_excel(writer, sheet_name="Ranking of Samples", index=False)
            for metric, pvals in pairwise_all.items():
                rows = []
                for (s1, s2), p in pvals.items():
                    rows.append({"Sample 1": s1, "Sample 2": s2, "p-value": p, "Excel Calculation": "=T.TEST(array1,array2,2,3)"})
                df_metric = pd.DataFrame(rows)
                sheet_name = sanitize_sheetname(metric)
                df_metric.to_excel(writer, sheet_name=sheet_name, index=False)
            for day in days_list:
                df_day = pd.DataFrame({
                    "Sample Casted on Electrode": aggregated_df["Sample Casted on Electrode"],
                    f"Sensitivity on Day {day} (mA/M)": aggregated_df[f"Sensitivity on Day {day} (mA/M)"],
                    f"Std Dev on Day {day} (mA/M)": std_df[f"Standard Deviation of Sensitivity on Day {day} (mA/M)"]
                })
                sheet_name = sanitize_sheetname(f"Day{day}_Bar")
                df_day.to_excel(writer, sheet_name=sheet_name, index=False)
            df_avg = pd.DataFrame({
                "Sample Casted on Electrode": aggregated_df["Sample Casted on Electrode"],
                "Average Sensitivity Across Days (mA/M)": aggregated_df["Average Sensitivity Across Days (mA/M)"],
                "Std Dev in Average Sensitivity (mA/M)": std_df["Standard Deviation in Average Sensitivity Across Days (mA/M)"]
            })
            df_avg.to_excel(writer, sheet_name="Avg_Bar", index=False)
            combined_clustered_df = aggregated_df.copy()
            for day in days_list:
                combined_clustered_df[f"Std Dev on Day {day} (mA/M)"] = std_df[f"Standard Deviation of Sensitivity on Day {day} (mA/M)"].values
            combined_clustered_df.to_excel(writer, sheet_name="Combined_Clustered_Bar", index=False)
            scatter_data = {"Day": days_list}
            for sample in aggregated_df["Sample Casted on Electrode"]:
                medians = []
                for day in days_list:
                    vals = agg_dict[sample][day]
                    med = np.nanmedian(vals) if len(vals)>0 else np.nan
                    medians.append(med)
                scatter_data[sample] = medians
            df_scatter = pd.DataFrame(scatter_data)
            df_scatter.to_excel(writer, sheet_name="Combined_Scatter", index=False)
            df_scatter.to_excel(writer, sheet_name="Combined_Scatter_Smooth", index=False)
        
        # --- Generate Aggregated Graphs ---
        day_bar = {}
        aggregated_day_data = {}
        aggregated_annotation_data = {}
        for day in days_list:
            plt.figure()
            samples = aggregated_df["Sample Casted on Electrode"]
            sens_vals = aggregated_df[f"Sensitivity on Day {day} (mA/M)"]
            std_vals = []
            for cast in samples:
                std_val = std_df.loc[std_df["Sample Casted on Electrode"]==cast, f"Standard Deviation of Sensitivity on Day {day} (mA/M)"]
                std_vals.append(std_val.values[0] if not std_val.empty else 0)
            bars = plt.bar(samples, sens_vals, yerr=std_vals, capsize=5, color='skyblue')
            plt.title(textwrap.fill(f"Graph of Sensitivity of Electrodes Casted with Samples on Day {day}", width=40))
            plt.xlabel("Sample")
            plt.ylabel("Sensitivity (mA/M)")
            plt.xticks(rotation=45, ha='right')
            for i in range(len(samples)-1):
                s1 = samples.iloc[i]
                s2 = samples.iloc[i+1]
                p_val = pairwise_all.get(f"Sensitivity on Day {day} (mA/M)", {}).get((s1, s2), None)
                if p_val is not None and p_val < 0.05:
                    annotate_significance(plt.gca(), i, i+1, sens_vals.iloc[i]+std_vals[i], 0.05, "*")
            plt.tight_layout()
            bar_path = os.path.join(app.config['STATIC_FOLDER'], f"agg_day_{day}.png")
            plt.savefig(bar_path, bbox_inches="tight")
            plt.close()
            day_bar[day] = os.path.basename(bar_path)
            aggregated_day_data[day] = pd.DataFrame({
                "Sample Casted on Electrode": samples,
                f"Sensitivity on Day {day} (mA/M)": sens_vals,
                f"Std Dev on Day {day} (mA/M)": std_vals
            }).to_html(classes='small-data', index=False)
            aggregated_annotation_data[day] = "<p style='font-size:10px;'>Significance annotations applied where p < 0.05</p>"
        
        plt.figure()
        samples = aggregated_df["Sample Casted on Electrode"]
        avg_sens = aggregated_df["Average Sensitivity Across Days (mA/M)"]
        std_avg = []
        for cast in samples:
            std_val = std_df.loc[std_df["Sample Casted on Electrode"]==cast, "Standard Deviation in Average Sensitivity Across Days (mA/M)"]
            std_avg.append(std_val.values[0] if not std_val.empty else 0)
        plt.bar(samples, avg_sens, yerr=std_avg, capsize=5, color='green')
        plt.title(textwrap.fill("Graph of Average Sensitivity Across the Days of Electrodes Casted with Sample", width=40))
        plt.xlabel("Sample")
        plt.ylabel("Average Sensitivity (mA/M)")
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        avg_bar_path = os.path.join(app.config['STATIC_FOLDER'], "avg_bar.png")
        plt.savefig(avg_bar_path, bbox_inches="tight")
        plt.close()
        
        plt.figure()
        unique_samples = aggregated_df["Sample Casted on Electrode"].tolist()
        n_samples = len(unique_samples)
        n_days = len(days_list)
        x = np.arange(n_samples)
        width = 0.8 / n_days
        for idx, day in enumerate(days_list):
            sens = aggregated_df[f"Sensitivity on Day {day} (mA/M)"]
            std_vals = []
            for cast in unique_samples:
                std_val = std_df.loc[std_df["Sample Casted on Electrode"]==cast, f"Standard Deviation of Sensitivity on Day {day} (mA/M)"]
                std_vals.append(std_val.values[0] if not std_val.empty else 0)
            plt.bar(x + idx*width, sens, width=width, yerr=std_vals, capsize=3, label=f"Day {day}")
        plt.xticks(x + width*(n_days-1)/2, unique_samples, rotation=45, ha='right')
        plt.title(textwrap.fill("Combined Clustered Bar Graph for All Samples", width=40))
        plt.xlabel("Sample")
        plt.ylabel("Sensitivity (mA/M)")
        plt.legend()
        plt.tight_layout()
        clustered_bar_path = os.path.join(app.config['STATIC_FOLDER'], "clustered_bar.png")
        plt.savefig(clustered_bar_path, bbox_inches="tight")
        plt.close()
        combined_clustered_df = aggregated_df.copy()
        for day in days_list:
            combined_clustered_df[f"Std Dev on Day {day} (mA/M)"] = std_df[f"Standard Deviation of Sensitivity on Day {day} (mA/M)"].values
        combined_clustered_table_html = combined_clustered_df.to_html(classes='small-data', index=False)
        
        sample_day_bar = {}
        sample_day_data = {}
        for sample in unique_samples:
            x_positions = list(range(len(days_list)))
            y_vals = []
            y_err = []
            table_rows = []
            for day in days_list:
                vals = agg_dict[sample][day]
                avg_val = np.nanmean(vals) if len(vals)>0 else np.nan
                std_val = np.std(vals, ddof=1) if len(vals)>1 else 0
                y_vals.append(avg_val)
                y_err.append(std_val)
                table_rows.append({"Day": day, "Sensitivity (mA/M)": avg_val, "Std Dev (mA/M)": std_val})
            plt.figure()
            plt.bar(x_positions, y_vals, yerr=y_err, capsize=5, color='purple')
            plt.xticks(x_positions, days_list)
            plt.title(textwrap.fill(f"Graph of Sensitivity of Electrodes Casted with {sample} Across the Days", width=40))
            plt.xlabel("Day")
            plt.ylabel("Sensitivity (mA/M)")
            plt.tight_layout()
            sample_bar_path = os.path.join(app.config['STATIC_FOLDER'], f"sample_day_bar_{sample}.png")
            plt.savefig(sample_bar_path, bbox_inches="tight")
            plt.close()
            sample_day_bar[sample] = os.path.basename(sample_bar_path)
            sample_day_data[sample] = pd.DataFrame(table_rows).to_html(classes='small-data', index=False)
        
        sample_scatter_images = {}
        sample_scatter_smooth_images = {}
        sample_scatter_tables = {}
        sample_scatter_smooth_tables = {}
        
        for sample in unique_samples:
            x_points = []
            y_points = []
            table_rows = []
            for day in days_list:
                try:
                    day_num = float(day)
                except:
                    day_num = day
                x_points.append(day_num)
                vals = agg_dict[sample][day]
                median_val = np.nanmedian(vals) if len(vals)>0 else np.nan
                y_points.append(median_val)
                table_rows.append({"Day": day, "Median Sensitivity (mA/M)": median_val})
            plt.figure()
            plt.scatter(x_points, y_points, label=sample)
            if len(x_points) > 1:
                coeffs = np.polyfit(x_points, y_points, 1)
                trendline = np.poly1d(coeffs)
                x_sorted = np.sort(x_points)
                y_fit = trendline(x_sorted)
                ss_res = np.sum((np.array(y_points) - trendline(np.array(x_points)))**2)
                ss_tot = np.sum((np.array(y_points) - np.mean(y_points))**2)
                r2_val = 1 - ss_res/ss_tot if ss_tot!=0 else 0
                plt.plot(x_sorted, y_fit, '--', label=f"Trend: y={coeffs[0]:.4f}x+{coeffs[1]:.4f}, R²={r2_val:.4f}")
            plt.title(textwrap.fill(f"Graph of Sensitivity of Electrodes Casted with {sample} Across the Days (Scatter Plot)", width=40))
            plt.xlabel("Day")
            plt.ylabel("Sensitivity (mA/M)")
            plt.legend(loc='best')
            plt.tight_layout()
            scatter_path = os.path.join(app.config['STATIC_FOLDER'], f"scatter_{sample}.png")
            plt.savefig(scatter_path, bbox_inches="tight")
            plt.close()
            sample_scatter_images[sample] = os.path.basename(scatter_path)
            sample_scatter_tables[sample] = pd.DataFrame(table_rows).to_html(classes='small-data', index=False)
            
            plt.figure()
            electrode_days = {}
            for elec in electrode_infos:
                if elec['cast'] == sample:
                    electrode_days[elec['actual']] = []
                    for day in days_list:
                        electrode_days[elec['actual']].append(sensitivity_dict[elec['actual']].get(day, np.nan))
            table_rows2 = []
            for elec_num, sens_list in electrode_days.items():
                try:
                    x_numeric = [float(day) for day in days_list]
                except:
                    x_numeric = list(range(len(days_list)))
                plt.scatter(x_numeric, sens_list, label=f"Electrode {elec_num}")
                if len(x_numeric)>1:
                    coeffs = np.polyfit(x_numeric, sens_list, 1)
                    trendline = np.poly1d(coeffs)
                    x_sorted = np.sort(x_numeric)
                    y_fit = trendline(x_sorted)
                    plt.plot(x_sorted, y_fit, '-', label=f"Smooth: Elec {elec_num}")
                    table_rows2.append({"Electrode": elec_num, "Gradient": coeffs[0]})
            plt.title(textwrap.fill(f"Graph of Sensitivity of Electrodes Casted with {sample} Across the Days (Smooth Curve)", width=40))
            plt.xlabel("Day")
            plt.ylabel("Sensitivity (mA/M)")
            plt.legend(loc='best')
            plt.tight_layout()
            scatter_smooth_path = os.path.join(app.config['STATIC_FOLDER'], f"scatter_smooth_{sample}.png")
            plt.savefig(scatter_smooth_path, bbox_inches="tight")
            plt.close()
            sample_scatter_smooth_images[sample] = os.path.basename(scatter_smooth_path)
            sample_scatter_smooth_tables[sample] = pd.DataFrame(table_rows2).to_html(classes='small-data', index=False)
        
        plt.figure()
        for sample in unique_samples:
            x_points = []
            y_points = []
            for day in days_list:
                try:
                    day_num = float(day)
                except:
                    day_num = day
                x_points.append(day_num)
                vals = agg_dict[sample][day]
                median_val = np.nanmedian(vals) if len(vals)>0 else np.nan
                y_points.append(median_val)
            plt.scatter(x_points, y_points, label=sample)
            if len(x_points)>1:
                coeffs = np.polyfit(x_points, y_points, 1)
                trendline = np.poly1d(coeffs)
                x_sorted = np.sort(x_points)
                y_fit = trendline(x_sorted)
                ss_res = np.sum((np.array(y_points)-trendline(np.array(x_points)))**2)
                ss_tot = np.sum((np.array(y_points)-np.mean(y_points))**2)
                r2_val = 1 - ss_res/ss_tot if ss_tot!=0 else 0
                plt.plot(x_sorted, y_fit, '--', label=f"{sample} Fit: R²={r2_val:.4f}")
        plt.title(textwrap.fill("Graph of Sensitivity of Electrodes Casted with Samples Across the Days (Combined Scatter Plot)", width=40))
        plt.xlabel("Day")
        plt.ylabel("Sensitivity (mA/M)")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        combined_scatter_path = os.path.join(app.config['STATIC_FOLDER'], "combined_scatter.png")
        plt.savefig(combined_scatter_path, bbox_inches="tight")
        plt.close()
        
        plt.figure()
        for sample in unique_samples:
            x_points = []
            y_points = []
            for day in days_list:
                try:
                    day_num = float(day)
                except:
                    day_num = day
                x_points.append(day_num)
                vals = agg_dict[sample][day]
                median_val = np.nanmedian(vals) if len(vals)>0 else np.nan
                y_points.append(median_val)
            plt.scatter(x_points, y_points, label=sample)
            if len(x_points) > 1:
                x_sorted = np.sort(x_points)
                y_sorted = []
                for xx in x_sorted:
                    vals_for_day = [v for (d,v) in zip(x_points,y_points) if d==xx]
                    y_sorted.append(np.nanmedian(vals_for_day))
                plt.plot(x_sorted, y_sorted, '-', label=f"{sample} Smooth")
                coeffs = np.polyfit(x_points, y_points, 1)
                trendline = np.poly1d(coeffs)
                plt.plot(x_sorted, trendline(x_sorted), ':', label=f"{sample} Dotted: y={coeffs[0]:.4f}x+{coeffs[1]:.4f}, R²={r2_val:.4f}")
        plt.title(textwrap.fill("Graph of Sensitivity of Electrodes Casted with Samples Across the Days (Combined Scatter with Smooth Curve)", width=40))
        plt.xlabel("Day")
        plt.ylabel("Sensitivity (mA/M)")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        combined_scatter_smooth_path = os.path.join(app.config['STATIC_FOLDER'], "combined_scatter_smooth.png")
        plt.savefig(combined_scatter_smooth_path, bbox_inches="tight")
        plt.close()
        
        # --- Generate Excel file for p-value details for all metrics ---
        pvalues_excel_path = os.path.join(app.config['STATIC_FOLDER'], "pvalues_summary.xlsx")
        with pd.ExcelWriter(pvalues_excel_path, engine='xlsxwriter') as writer:
            for metric, pvals in pairwise_all.items():
                rows = []
                for (s1, s2), p in pvals.items():
                    rows.append({"Sample 1": s1, "Sample 2": s2, "p-value": p, "Excel Formula": "=T.TEST(array1,array2,2,3)"})
                df_metric = pd.DataFrame(rows)
                sheet_name = sanitize_sheetname(metric)
                df_metric.to_excel(writer, sheet_name=sheet_name, index=False)
            explanation = ("Pairwise Welch's t-tests (unequal variances) were performed for each sensitivity metric. "
                           "Samples with p-value < 0.05 are considered significantly different; otherwise, they share the same rank.")
            ws = writer.book.add_worksheet("Explanation")
            ws.write(0, 0, explanation)
        
        # --- Consolidate Individual Electrode Excel Files with Raw and Processed Data ---
        combined_indiv_excel_path = os.path.join(app.config['STATIC_FOLDER'], "Summary_of_Individual_Electrode_Electrochemical_Sensitivity.xlsx")
        with pd.ExcelWriter(combined_indiv_excel_path, engine='xlsxwriter') as writer:
            for key, result in individual_results.items():
                sheet_name_raw = sanitize_sheetname(f"Elec{result['electrode']}_Day{result['day']}_Raw")
                df_raw = pd.read_html(result['table_html'])[0]
                df_raw.to_excel(writer, sheet_name=sheet_name_raw, index=False)
                worksheet_raw = writer.sheets[sheet_name_raw]
                raw_graph_path_full = os.path.join(app.config['STATIC_FOLDER'], result['raw_graph'])
                worksheet_raw.insert_image('G2', raw_graph_path_full)
                raw_graph_500_path_full = os.path.join(app.config['STATIC_FOLDER'], result['raw_graph_500'])
                worksheet_raw.insert_image('G20', raw_graph_500_path_full)
                
                sheet_name_proc = sanitize_sheetname(f"Elec{result['electrode']}_Day{result['day']}_Proc")
                df_proc = pd.read_html(result['table_html'])[0]
                df_proc.to_excel(writer, sheet_name=sheet_name_proc, index=False)
                worksheet_proc = writer.sheets[sheet_name_proc]
                processed_graph_path_full = os.path.join(app.config['STATIC_FOLDER'], result['processed_graph'])
                worksheet_proc.insert_image('G2', processed_graph_path_full)
        
        scatter_data = {"Day": days_list}
        for sample in aggregated_df["Sample Casted on Electrode"]:
            medians = []
            for day in days_list:
                vals = agg_dict[sample][day]
                med = np.nanmedian(vals) if len(vals)>0 else np.nan
                medians.append(med)
            scatter_data[sample] = medians
        df_scatter = pd.DataFrame(scatter_data)
        combined_scatter_table = df_scatter.to_html(classes='small-data', index=False)
        
        downloads = {
            "combined_excel": os.path.basename(combined_excel_path),
            "pvalues_excel": os.path.basename(pvalues_excel_path),
            "avg_bar": os.path.basename(avg_bar_path),
            "clustered_bar": os.path.basename(clustered_bar_path),
            "combined_scatter": os.path.basename(combined_scatter_path),
            "combined_scatter_smooth": os.path.basename(combined_scatter_smooth_path),
            "combined_indiv_excel": os.path.basename(combined_indiv_excel_path)
        }
        
        return render_template(
                'result.html',
                days_list=days_list,
                individual_table_html=summary_df.to_html(classes='data', border=0, index=False),
                aggregated_table_html=aggregated_df.to_html(classes='data', border=0, index=False),
                std_table_html=std_df.to_html(classes='data', border=0, index=False),
                ranking_table=ranking_table_html,
                pvalue_tables_html=pvalues_tables_html,
                combined_clustered_table=combined_clustered_table_html,
                day_bar=day_bar,
                aggregated_day_data=aggregated_day_data,
                aggregated_annotation_data=aggregated_annotation_data,
                avg_bar=os.path.basename(avg_bar_path),
                clustered_bar=os.path.basename(clustered_bar_path),
                sample_day_bar=sample_day_bar,
                sample_day_data=sample_day_data,
                sample_scatter=sample_scatter_images,
                sample_scatter_tables=sample_scatter_tables,
                sample_scatter_smooth=sample_scatter_smooth_images,
                sample_scatter_smooth_tables=sample_scatter_smooth_tables,
                combined_scatter=os.path.basename(combined_scatter_path),
                combined_scatter_smooth=os.path.basename(combined_scatter_smooth_path),
                combined_scatter_table=combined_scatter_table,
                downloads=downloads,
                individual_results=individual_results
            )

    return render_template('index.html')

@app.route('/download/<filename>')
def download_file(filename):
    return send_from_directory(app.config['STATIC_FOLDER'], filename, as_attachment=True)

if __name__ == '__main__':
    app.run(debug=True)


# In[ ]:




