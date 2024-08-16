import os
from statsmodels.stats.multitest import multipletests
from scipy.stats import f_oneway as anova
import scipy
from scipy.stats import t
import plotly.express as px
from pathlib import Path
  
# from select import select
import numpy as np
import pandas as pd

import plotly.graph_objs as go
# from plotly.offline import plot
# import plotly.io as pio
# from plotly.offline import iplot
# import plotly as py
from matplotlib_venn import venn2, venn3
import matplotlib.pyplot as plt
# import matplotlib
# import plotly.colors
import plotly.express as px
# from plotly.subplots import make_subplots

from SCP_processor import SCP_processor


class SCP_plotter:
    def __init__(self, write_output = False, data_type = "LF", processor = SCP_processor()) -> None:
        self.write_output = write_output
        self.app_folder = "./"
        self.processor = processor
        self.url_base = None
        self.data_type = data_type


    def ID_plots(self, data_object, plot_options, saved_settings, username=None):
        """_Prepare data for creating protein peptide identification bar
        plot_

        Args:
            data_dict (_type_): _description_
        """
        # Create an empty dictionary to store the group names and filters
        group_names = [key for key in saved_settings.keys() if "Order@" not in str(key)]

        # import the data and save order
        group_dict = {}

        if plot_options["ID mode"] == "MS2" or plot_options["ID mode"] == "total" or plot_options["ID mode"] == "stacked":
            x_axis_order = saved_settings["Order@Conditions"]
        elif plot_options["Group By X"] == "ID_Mode":
            print("ERROR: x axis separation of MS2/MBR not supported")
        else:
            x_axis_order = saved_settings["Order@"+plot_options["Group By X"]]
        if plot_options["ID mode"] == "grouped" or plot_options["ID mode"] == "grouped_stacked" and plot_options["Group By Color"] != "ID_Mode":
            color_order = saved_settings["Order@"+plot_options["Group By Color"]]
            plot_options["color_order"] = color_order

        plot_options["x_axis_order"] = x_axis_order

        # filter runs into different groups
        i = 1
        runname_list = []  # contain list of run names list for each groups
        for eachGroup in group_names:
            run_ids = saved_settings[eachGroup]["records"]
            runname_sublist = data_object["run_metadata"].loc[data_object["run_metadata"]["Run Identifier"].isin(run_ids)]["Run Names"].to_list()
            group_dict[eachGroup] =  self.processor.filter_by_id(
                    data_object,
                    run_ids) # prevent the list from being changed
            runname_list.append(runname_sublist)
            i += 1

            #print(group_dict[eachGroup]["run_metadata"])
        #display(data_object["protein_ID_Summary"])
        #display(group_dict[eachGroup]["protein_ID_Summary"])
        # create ID plots
        # allIDs table will be used to store all experiment name, ID types (
        # protein, peptide, MS2 and MS1 based), conditions and IDs numbers
        allIDs = pd.DataFrame(
            columns=["Names", "ID_Type", "ID_Mode", "Conditions", "IDs"])

        # loop through each group and extract IDs, put them into allIDs table
        
            
        # ######################allIDs format###################
        # name	ID_Type	ID_Mode	Conditions	IDs
        # file1	peptide	MS2_IDs	experimetn 1	xxxxx
        # file2	protein	MBR_IDs	experiment 2	xxxx
        # file3	peptide	Total_IDs	experiment 3	xxx
        #######################################################
        # Calcuate mean, standard deviation and number of replicates for each
        # choose protein or peptide
        if plot_options["plot_type"] == "1":  # Protein ID
            for eachCondition in group_names:
            # Protein ID summary
                for index, row in group_dict[eachCondition][
                        "protein_ID_Summary"].iterrows():
                    for item in ["MS2_IDs",
                                "MBR_IDs",
                                "Total_IDs"]:
                        if not pd.isna(group_dict[eachCondition][
                                "protein_ID_Summary"].at[index, item]):
                            # if the row with the item column is not empty,
                            # add it to allIDs table.
                            allIDs = pd.concat(
                                [allIDs,
                                pd.DataFrame(
                                    [[group_dict[eachCondition][
                                        "protein_ID_Summary"].at[index, "names"],
                                    "protein",
                                    item,
                                    eachCondition,
                                    group_dict[eachCondition][
                                        "protein_ID_Summary"].at[index, item]]],
                                    columns=["Names",
                                            "ID_Type",
                                            "ID_Mode",
                                            "Conditions",
                                            "IDs"])],
                                ignore_index=True)
            allIDs = allIDs[allIDs["ID_Type"] == "protein"]
                            
        elif plot_options["plot_type"] == "2":  # Peptide ID
            # Peptide ID summary
            for eachCondition in group_names:

                for index, row in group_dict[eachCondition][
                        "peptide_ID_Summary"].iterrows():
                    for item in ["MS2_IDs",
                                "MBR_IDs",
                                "Total_IDs"]:
                        if not pd.isna(group_dict[eachCondition][
                                "peptide_ID_Summary"].at[index, item]):
                            allIDs = pd.concat(
                                [allIDs,
                                pd.DataFrame(
                                    [[group_dict[eachCondition][
                                        "peptide_ID_Summary"].at[index, "names"],
                                    "peptide",
                                    item,
                                    eachCondition,
                                    group_dict[eachCondition][
                                        "peptide_ID_Summary"].at[index, item]]],
                                    columns=["Names",
                                            "ID_Type",
                                            "ID_Mode",
                                            "Conditions",
                                            "IDs"])],
                                ignore_index=True)
            allIDs = allIDs[allIDs["ID_Type"] == "peptide"]
                        
        export_ids = allIDs.copy()

        # choose total, MS2 or stacked
        if plot_options["ID mode"] == "MS2":  # MS2 ID
            allIDs = allIDs[allIDs["ID_Mode"] == "MS2_IDs"]
        elif plot_options["ID mode"] == "total":
            # total ID combined, if not already summed (key exist), sum them
    #         if allIDs[allIDs["ID_Mode"] == "Total_IDs"].empty:
    #             grouped = allIDs.groupby('name').agg(
    #                 {'IDs': 'sum', 'ID_Type': 'first', 'Conditions': 'first'})
    #             grouped = grouped.reset_index()
    #             grouped["ID_Mode"] = "Total_IDs"
    #             allIDs = grouped
            allIDs = allIDs[allIDs["ID_Mode"] == "Total_IDs"]
        elif plot_options["Group By X"] == "ID_Mode" or plot_options["Group By Color"] == "ID_Mode" \
            or plot_options["Group By Stack"] == "ID_Mode" and not (plot_options["ID mode"] == "total" or plot_options["ID mode"] == "MS2"):  # total separated
            pass
        elif plot_options["ID mode"] == "MS2":
            allIDs = allIDs[allIDs["ID_Mode"] == "MS2_IDs"]
        else:
            allIDs = allIDs[allIDs["ID_Mode"] == "Total_IDs"]

        toPlotIDs = allIDs.groupby(["ID_Mode", "Conditions"]).agg({
            'IDs': ['mean', 'std', 'count'], 'ID_Type': 'first', })

        # rename the columns
        toPlotIDs.columns = ['IDs', 'stdev', 'n', 'ID_Type']
        # reset the index after grouping
        toPlotIDs = toPlotIDs.reset_index()
        # calculate the confidence interval based on 95%confidence interval`
        toPlotIDs["confInt"] = t.ppf(0.975, toPlotIDs['n']-1) * \
            toPlotIDs['stdev']/np.sqrt(toPlotIDs['n'])

        #add columns for the categories specified in settings file (the one with all the filenames)
        standard_groups = ["filter_in","filter_out","records"]
        categories = [col for col in list(saved_settings[list(group_names)[0]].keys()) if col not in standard_groups]
        for eachCategory in categories:
            toPlotIDs[eachCategory] = ""
            for eachGroup in group_names:
                toPlotIDs.loc[toPlotIDs["Conditions"]==eachGroup,eachCategory] = saved_settings[eachGroup][eachCategory]
        
        #display(toPlotIDs)
        fig = self.plot_IDChart_plotly(toPlotIDs,
                                username=username,
                                plot_options=plot_options)

        if self.write_output:    
            # export the data to csv for user downloading
            data_dir = os.path.join(self.app_folder, "csv/")
            # create the directory if it does not exist
            if not os.path.exists(data_dir):
                Path(data_dir).mkdir(parents=True)
            categories = [col for col in list(saved_settings[list(group_names)[0]].keys()) if col not in standard_groups]
            
            for eachCategory in categories:
                export_ids[eachCategory] = ""
                for eachGroup in group_names:
                    export_ids.loc[export_ids["Conditions"]==eachGroup,eachCategory] = saved_settings[eachGroup][eachCategory]
            export_ids = export_ids.replace({"Names": dict(zip(data_object["run_metadata"]["Run Identifier"],data_object["run_metadata"]["Run Names"]))})
            
            # export the data to csv
            export_ids.to_csv(os.path.join(
                data_dir, f"{username}_ID_data.csv"), index=False)
            
            print("Downloading links...")

            # create the link for downloading the data
            CSV_link = f"/files/{self.url_base}/csv/" \
                f"{username}_ID_data.csv"

            # add SVG download link

            SVG_link = f"/files/{self.url_base}/images/" \
                f"{username}_ID_Bar_Plot.svg"

            img_dir = os.path.join(self.app_folder, "images/")
            if not os.path.exists(img_dir):
                Path(img_dir).mkdir(parents=True)

            fig.write_image(os.path.join(
                img_dir, f"{username}_ID_Bar_Plot.svg"), format = "svg", validate = False, engine = "kaleido")


        else:
            CSV_link = None
            SVG_link = None
        return fig, CSV_link, SVG_link

    def plot_IDChart_plotly(self, ID_data,
                            username=None,
                            plot_options=None):
        """_Plot the ID bar plot for the given data_

        Args:
            ID_data (_type_): _description_
            username (str, optional): _description_. Defaults to "test".
            plot_options (_type_, optional): _description_. Defaults to None.

        Returns:
            _type_: _description_
        """

        plot_div = None
        

        if plot_options["ID mode"] == "grouped":  
            # plot options
            # error bar
            if plot_options["error bar"] == "stdev":
                error_bars = "stdev"
                error_visibile = True
            elif plot_options["error bar"] == "ci95":
                error_bars = "confInt"
                error_visibile = True
            else:
                error_bars = "stdev"
                error_visibile = False

            # mean label
            if plot_options["mean label"] == "True" or \
                    plot_options["mean label"] == True:
                total_labels = [{"x": x, "y": total*1.15, "text": str(
                    int(total)), "showarrow": False} for x, total in zip(
                        ID_data["Conditions"], ID_data["IDs"])]
            else:
                total_labels = []   # no mean labels

            if plot_options["Group By X"] == "ID_Mode" or plot_options["Group By Color"] == "ID_Mode":  # total separated
                ID_data = ID_data[ID_data["ID_Mode"] != "Total_IDs"]
            else:
                ID_data = ID_data[ID_data["ID_Mode"] == "Total_IDs"]
            #find out present categories
            categories = plot_options["color_order"]
            # create the plot
            fig_data = []
            i = 0
            for eachCategory in categories:
                fig_data.append(go.Bar(name = eachCategory,
                            x=ID_data.loc[ID_data[plot_options["Group By Color"]]==eachCategory,plot_options["Group By X"]].tolist(),
                            y=ID_data.loc[ID_data[plot_options["Group By Color"]]==eachCategory,"IDs"].tolist(),
                            marker_color = plot_options["color"][i],
                            text = [round(x) for x in ID_data.loc[ID_data[plot_options["Group By Color"]]==eachCategory,"IDs"].tolist()],
                            error_y=dict(
                                type = "data",
                                array = ID_data.loc[ID_data[plot_options["Group By Color"]]==eachCategory,error_bars].tolist(),
                                visible = error_visibile
                            )))
                i = i + 1                    

            fig = go.Figure(data = fig_data,
                            layout=go.Layout(yaxis_title=plot_options["Y Title"],
                            xaxis_title=plot_options["Group By X"],
                            barmode="group",paper_bgcolor="rgba(255,255,255,255)",
                            plot_bgcolor="rgba(255, 255, 255, 255)",
                            yaxis=dict(showline=True, linewidth=1, linecolor='black'),
                            xaxis=dict(showline=True, linewidth=1, linecolor='black')))
            fig.update_xaxes(categoryorder='array', categoryarray = plot_options["x_axis_order"])

        elif plot_options["ID mode"] == "stacked":
            # plot options
            # error bar
            if plot_options["error bar"] == "stdev":
                error_bars = "stdev"
                error_visibile = True
            elif plot_options["error bar"] == "ci95":
                error_bars = "confInt"
                error_visibile = True
            else:
                error_bars = "stdev"
                error_visibile = False

            # mean label
            if plot_options["mean label"] == "True" or \
                    plot_options["mean label"] == True:
                total_labels = [{"x": x, "y": total*1.15, "text": str(
                    int(total)), "showarrow": False} for x, total in zip(
                        ID_data["Conditions"], ID_data["IDs"])]
            else:
                total_labels = []   # no mean labels
            if plot_options["Group By X"] == "ID_Mode" or plot_options["Group By Stack"] == "ID_Mode":  # total separated
                ID_data = ID_data[ID_data["ID_Mode"] != "Total_IDs"]
            else:
                ID_data = ID_data[ID_data["ID_Mode"] == "Total_IDs"]

            if plot_options["Group By Stack"] == "ID_Mode":
                layers = ["MS2_IDs", "MBR_IDs"]
            else:
                layers = ID_data.groupby(plot_options["Group By Stack"]).first().reset_index()[plot_options["Group By Stack"]].tolist()
            fig_data = []
            last_layer = None
            i = 0
            for eachLayer in layers: 
                if last_layer == None:
                    fig_data.append(go.Bar(
                        name = eachLayer,
                        x = ID_data.loc[(ID_data[plot_options["Group By Stack"]]==eachLayer),plot_options["Group By X"]].tolist(),
                        y = ID_data.loc[(ID_data[plot_options["Group By Stack"]]==eachLayer),"IDs"].tolist(),
                        marker_color = plot_options["color"][i],
                        text = [round(x) for x in ID_data.loc[(ID_data[plot_options["Group By Stack"]]==eachLayer),"IDs"].tolist()],
                        error_y= dict(
                            type = "data",
                            array = ID_data.loc[ID_data[plot_options["Group By Stack"]]==eachLayer,error_bars].tolist(),
                            visible = error_visibile
                        )
                    ))
                    bases = ID_data.loc[(ID_data[plot_options["Group By Stack"]]==eachLayer),"IDs"]
                    i = i + 1
                    
                else:
                    fig_data.append(go.Bar(
                        name = eachLayer,
                        x = ID_data.loc[(ID_data[plot_options["Group By Stack"]]==eachLayer),plot_options["Group By X"]].tolist(),
                        y = ID_data.loc[(ID_data[plot_options["Group By Stack"]]==eachLayer),"IDs"].tolist(),
                        base=bases,
                        marker_color = plot_options["color"][i],
                        opacity=0.5,
                        text = [round(x) for x in bases + ID_data.loc[(ID_data[plot_options["Group By Stack"]]==eachLayer),"IDs"].tolist()],
                        error_y= dict(
                            type = "data",
                            array = ID_data.loc[ID_data[plot_options["Group By Stack"]]==eachLayer,error_bars].tolist(),
                            visible = error_visibile
                        )
                    ))
                    print(bases)
                    bases = bases + ID_data.loc[(ID_data[plot_options["Group By Stack"]]==eachLayer),"IDs"].tolist()
                last_layer = eachLayer
            fig = go.Figure(
                    data = fig_data,
                    layout=go.Layout(
                    yaxis_title=plot_options["Y Title"],
                    xaxis_title=plot_options["Group By X"],
                    barmode="stack", 
                    paper_bgcolor="rgba(255,255,255,255)",
                    plot_bgcolor="rgba(255, 255, 255, 255)",
                    yaxis=dict(showline=True, linewidth=1, linecolor='black'),
                    xaxis=dict(showline=True, linewidth=1, linecolor='black')
                ))          
            fig.update_xaxes(categoryorder='array', categoryarray = plot_options["x_axis_order"])            
        elif plot_options["ID mode"] == "grouped_stacked":
            # plot options
            # error bar
            if plot_options["error bar"] == "stdev":
                error_bars = "stdev"
                error_visibile = True
            elif plot_options["error bar"] == "ci95":
                error_bars = "confInt"
                error_visibile = True
            else:
                error_bars = "stdev"
                error_visibile = False

            # mean label
            if plot_options["mean label"] == "True" or \
                    plot_options["mean label"] == True:
                total_labels = [{"x": x, "y": total*1.15, "text": str(
                    int(total)), "showarrow": False} for x, total in zip(
                        ID_data["Conditions"], ID_data["IDs"])]
            else:
                total_labels = []   # no mean labels

            if plot_options["Group By X"] == "ID_Mode" or plot_options["Group By Color"] == "ID_Mode"or plot_options["Group By Stack"] == "ID_Mode":  # total separated
                ID_data = ID_data[ID_data["ID_Mode"] != "Total_IDs"]
            else:
                ID_data = ID_data[ID_data["ID_Mode"] == "Total_IDs"]

            #make data tidy
            if plot_options["Group By Stack"] == "ID_Mode":
                layers = ["MS2_IDs", "MBR_IDs"]
            else:
                layers = ID_data.groupby(plot_options["Group By Stack"]).first().reset_index()[plot_options["Group By Stack"]].tolist()
            categories = plot_options["color_order"]
            
            fig_data = []
            i = 0
            for eachCategory in categories:
                last_layer = None
                j = 0
                for eachLayer in layers: 
                    if last_layer == None:
                        fig_data.append(go.Bar(
                            name = str(eachLayer) + " " + str(eachCategory),
                            x = ID_data.loc[(ID_data[plot_options["Group By Color"]]==eachCategory)&(ID_data[plot_options["Group By Stack"]]==eachLayer),plot_options["Group By X"]],
                            y = ID_data.loc[(ID_data[plot_options["Group By Color"]]==eachCategory)&(ID_data[plot_options["Group By Stack"]]==eachLayer),"IDs"],
                            offsetgroup=i,
                            text = [round(x) for x in ID_data.loc[(ID_data[plot_options["Group By Color"]]==eachCategory)&(ID_data[plot_options["Group By Stack"]]==eachLayer),"IDs"].tolist()],
                            marker_color = plot_options["color"][i],
                            error_y = dict(
                                type = "data",
                                array = ID_data.loc[(ID_data[plot_options["Group By Color"]]==eachCategory)&(ID_data[plot_options["Group By Stack"]]==eachLayer),error_bars],
                                visible=True)
                        ))
                        bases = ID_data.loc[(ID_data[plot_options["Group By Color"]]==eachCategory)&(ID_data[plot_options["Group By Stack"]]==eachLayer),"IDs"]
                    else:

                        fig_data.append(go.Bar(
                            name = str(eachLayer) + " " + str(eachCategory),
                            x = ID_data.loc[(ID_data[plot_options["Group By Color"]]==eachCategory)&(ID_data[plot_options["Group By Stack"]]==eachLayer),plot_options["Group By X"]],
                            y = ID_data.loc[(ID_data[plot_options["Group By Color"]]==eachCategory)&(ID_data[plot_options["Group By Stack"]]==eachLayer),"IDs"],
                            base=bases,
                            offsetgroup=i,
                            text = [round(x) for x in ID_data.loc[(ID_data[plot_options["Group By Color"]]==eachCategory)&(ID_data[plot_options["Group By Stack"]]==eachLayer),"IDs"].tolist()+bases],
                            marker_color = plot_options["color"][i],
                            opacity=1/2**j,
                            error_y = dict(
                                type = "data",
                                array = ID_data.loc[(ID_data[plot_options["Group By Color"]]==eachCategory)&(ID_data[plot_options["Group By Stack"]]==eachLayer),error_bars],
                                visible=True)
                            ))
                        bases = bases + ID_data.loc[(ID_data[plot_options["Group By Color"]]==eachCategory)&(ID_data[plot_options["Group By Stack"]]==eachLayer),"IDs"].tolist()
                    last_layer=eachLayer
                    j = j + 1
                i = i + 1

            fig = go.Figure(
                fig_data,
                layout=go.Layout(
                    yaxis_title=plot_options["Y Title"],
                    xaxis_title=plot_options["Group By X"],
                    barmode="group",
                    plot_bgcolor="rgba(255, 255, 255, 255)",
                    paper_bgcolor="rgba(255, 255, 255, 255)",
                    yaxis=dict(showline=True, linewidth=1, linecolor='black'),
                    xaxis=dict(showline=True, linewidth=1, linecolor='black')
                )
            )
            fig.update_xaxes(categoryorder='array', categoryarray = plot_options["x_axis_order"])
        else:

            # plot options
            # error bar
            if plot_options["error bar"] == "stdev":
                error_bars = "stdev"
            elif plot_options["error bar"] == "ci95":
                error_bars = "confInt"
            else:
                error_bars = None

            # mean label
            if plot_options["mean label"] == "True" or \
                    plot_options["mean label"] == True:
                total_labels = [{"x": x, "y": total*1.15, "text": str(
                    int(total)), "showarrow": False} for x, total in zip(
                        ID_data["Conditions"], ID_data["IDs"])]
            else:
                total_labels = []   # no mean labels

            # create the plot
            fig = px.bar(ID_data,
                        x="Conditions",
                        y="IDs",
                        error_y=error_bars,
                        color="Conditions",
                        color_discrete_sequence=plot_options["color"],
                        width=plot_options["width"],
                        height=plot_options["height"],
                        )
            fig.update_layout(xaxis_title=plot_options["X Title"],
                            yaxis_title=plot_options["Y Title"],
                            annotations=total_labels,
                            font=plot_options["font"],
                            plot_bgcolor ="rgba(255, 255, 255, 255)",
                            paper_bgcolor="rgba(255, 255, 255, 255)",
                            yaxis=dict(showline=True, linewidth=1, linecolor='black'),
                            xaxis=dict(showline=True, linewidth=1, linecolor='black')
                            )
            fig.update_xaxes(categoryorder='array', categoryarray = plot_options["x_axis_order"])
        return fig

    def Missing_Values_Plots(self, data_object, plot_options, saved_settings, username=None):
        """_Prepare data for creating protein peptide identification bar
        plot_

        Args:
            data_dict (_type_): _description_
        """
        # Create an empty dictionary to store the group names and filters
        group_names = [key for key in saved_settings.keys() if "Order@" not in str(key)]

        # import the data and save order
        group_dict = {}

        if plot_options["ID mode"] == "MS2" or plot_options["ID mode"] == "total" or plot_options["ID mode"] == "stacked":
            x_axis_order = saved_settings["Order@Conditions"]
        elif plot_options["Group By X"] == "ID_Mode":
            print("ERROR: x axis separation of MS2/MBR not supported")
        else:
            x_axis_order = saved_settings["Order@"+plot_options["Group By X"]]
        if plot_options["ID mode"] == "grouped" or plot_options["ID mode"] == "grouped_stacked" and plot_options["Group By Color"] != "ID_Mode":
            color_order = saved_settings["Order@"+plot_options["Group By Color"]]
            plot_options["color_order"] = color_order

        plot_options["x_axis_order"] = x_axis_order

        # filter runs into different groups
        i = 1
        runname_list = []  # contain list of run names list for each groups
        for eachGroup in group_names:
            run_ids = saved_settings[eachGroup]["records"]
            runname_sublist = data_object["run_metadata"].loc[data_object["run_metadata"]["Run Identifier"].isin(run_ids)]["Run Names"].to_list()
            if self.processor.data_type == "TMT":
                channel_ids = data_object["run_metadata"].loc[data_object["run_metadata"]["Run Identifier"].isin(run_ids)]["Channel Identifier"].to_list()
                group_dict[eachGroup] =  self.processor.filter_by_channel_id(
                    data_object,
                    channel_ids)  
            elif self.processor.data_type == "LF":

                group_dict[eachGroup] =  self.processor.filter_by_id(
                    data_object,
                    run_ids)  # prevent the list from being changed
            runname_list.append(runname_sublist)
            # print(runname_sublist)
            i += 1

            #print(group_dict[eachGroup]["run_metadata"])
        #display(data_object["protein_ID_Summary"])
        #display(group_dict[eachGroup]["protein_ID_Summary"])
        # create ID plots
        # allIDs table will be used to store all experiment name, ID types (
        # protein, peptide, MS2 and MS1 based), conditions and IDs numbers
        allProteins = pd.DataFrame(
            columns=["Accession", "Missing Values Rate","Conditions"])
        allPeptides = pd.DataFrame(
            columns=["Annotated Sequence", "Missing Values Rate","Conditions"])
        # loop through each group and extract IDs, put them into allIDs table
        for eachCondition in group_names:
            # Protein ID summary
            currentData = group_dict[eachCondition]
            if self.processor.ignore_proteins == False:
                current =  self.processor.calculate_missing_values_MS2(currentData, is_protein=True)
                # print(current)
                current["Conditions"] = eachCondition
                allProteins = pd.concat([allProteins, current])
                # print(allProteins)
            if self.processor.ignore_peptides == False:
                # Peptide ID summary
                current =  self.processor.calculate_missing_values_MS2(currentData, is_protein=False)
                current["Conditions"] = eachCondition
                allPeptides = pd.concat([allPeptides, current])

        
        # ######################allIDs format###################
        # name	ID_Type	ID_Mode	Conditions	IDs
        # file1	peptide	MS2_IDs	experimetn 1	xxxxx
        # file2	protein	MBR_IDs	experiment 2	xxxx
        # file3	peptide	Total_IDs	experiment 3	xxx
        #######################################################
        # Calcuate mean, standard deviation and number of replicates for each

        # choose protein or peptide
        if plot_options["plot_type"] == "1":  # Protein ID
            allIDs = allProteins
            name = "Accession"
        elif plot_options["plot_type"] == "2":  # Peptide ID
            allIDs = allPeptides
            name = "Annotated Sequence"


        toPlotIDs = pd.DataFrame(columns=["Conditions","Cutoff","IDs"])

        cutoffs = plot_options["cutoffs"]

        for eachCondition in group_names:
            last_cutoff = 0 #anything with more than 0 IDs
            for eachcutoff in cutoffs:
                currentIDs = allIDs.loc[(allIDs["Conditions"]==eachCondition) & (allIDs["Missing Values Rate"] >= last_cutoff) & (allIDs["Missing Values Rate"] < eachcutoff)].shape[0]
                # print(currentIDs)
                current = {"Conditions":eachCondition, "Cutoff": eachcutoff, "IDs": currentIDs}
                toPlotIDs = toPlotIDs.append(current, ignore_index = True)
                last_cutoff = eachcutoff

        
        #add columns for the categories specified in settings file (the one with all the filenames)
        standard_groups = ["filter_in","filter_out","records"]
        categories = [col for col in list(saved_settings[list(group_names)[0]].keys()) if col not in standard_groups]
        for eachCategory in categories:
            toPlotIDs[eachCategory] = ""
            for eachGroup in group_names:
                toPlotIDs.loc[toPlotIDs["Conditions"]==eachGroup,eachCategory] = saved_settings[eachGroup][eachCategory]
        
        #display(toPlotIDs)
        fig = self.plot_MissVal_plotly(toPlotIDs,
                                username=username,
                                plot_options=plot_options)


        if self.write_output:    
            # export the data to csv for user downloading
            data_dir = os.path.join(self.app_folder, "csv/")
            # create the directory if it does not exist
            if not os.path.exists(data_dir):
                Path(data_dir).mkdir(parents=True)
            categories = [col for col in list(saved_settings[list(group_names)[0]].keys()) if col not in standard_groups]
            
            
            toPlotIDs = toPlotIDs.replace({"Names": dict(zip(data_object["run_metadata"]["Run Identifier"],data_object["run_metadata"]["Run Names"]))})
            
            # export the data to csv
            toPlotIDs.to_csv(os.path.join(
                data_dir, f"{username}_MissVal_data.csv"), index=False)
            
            print("Downloading links...")

            # create the link for downloading the data
            CSV_link = f"/files/{self.url_base}/csv/" \
                f"{username}_MissVal_data.csv"

            # add SVG download link

            SVG_link = f"/files/{self.url_base}/images/" \
                f"{username}_MissVal_Bar_Plot.svg"

            img_dir = os.path.join(self.app_folder, "images/")
            if not os.path.exists(img_dir):
                Path(img_dir).mkdir(parents=True)

            fig.write_image(os.path.join(
                img_dir, f"{username}_MissVal_Bar_Plot.svg"), format = "svg", validate = False, engine = "kaleido")


        else:
            CSV_link = None
            SVG_link = None
        return fig, CSV_link, SVG_link

    def plot_MissVal_plotly(self, ID_data,
                            username=None,
                            plot_options=None):
        """_Plot the ID bar plot for the given data_

        Args:
            ID_data (_type_): _description_
            username (str, optional): _description_. Defaults to "test".
            plot_options (_type_, optional): _description_. Defaults to None.

        Returns:
            _type_: _description_
        """

        plot_div = None
        

        # if plot_options["ID mode"] == "grouped":  
        #     # plot options
        #     # error bar
        #     if plot_options["error bar"] == "stdev":
        #         error_bars = "stdev"
        #         error_visibile = True
        #     elif plot_options["error bar"] == "ci95":
        #         error_bars = "confInt"
        #         error_visibile = True
        #     else:
        #         error_bars = "stdev"
        #         error_visibile = False

        #     # mean label
        #     if plot_options["mean label"] == "True" or \
        #             plot_options["mean label"] == True:
        #         total_labels = [{"x": x, "y": total*1.15, "text": str(
        #             int(total)), "showarrow": False} for x, total in zip(
        #                 ID_data["Conditions"], ID_data["IDs"])]
        #     else:
        #         total_labels = []   # no mean labels

        #     if plot_options["Group By X"] == "ID_Mode" or plot_options["Group By Color"] == "ID_Mode":  # total separated
        #         ID_data = ID_data[ID_data["ID_Mode"] != "Total_IDs"]
        #     else:
        #         ID_data = ID_data[ID_data["ID_Mode"] == "Total_IDs"]
        #     #find out present categories
        #     categories = plot_options["color_order"]
        #     # create the plot
        #     fig_data = []
        #     i = 0
        #     for eachCategory in categories:
        #         fig_data.append(go.Bar(name = eachCategory,
        #                     x=ID_data.loc[ID_data[plot_options["Group By Color"]]==eachCategory,plot_options["Group By X"]].tolist(),
        #                     y=ID_data.loc[ID_data[plot_options["Group By Color"]]==eachCategory,"IDs"].tolist(),
        #                     marker_color = plot_options["color"][i],
        #                     text = [round(x) for x in ID_data.loc[ID_data[plot_options["Group By Color"]]==eachCategory,"IDs"].tolist()],
        #                     error_y=dict(
        #                         type = "data",
        #                         array = ID_data.loc[ID_data[plot_options["Group By Color"]]==eachCategory,error_bars].tolist(),
        #                         visible = error_visibile
        #                     )))
        #         i = i + 1                    

        #     fig = go.Figure(data = fig_data,
        #                     layout=go.Layout(yaxis_title=plot_options["Y Title"],
        #                     xaxis_title=plot_options["Group By X"],
        #                     barmode="group",paper_bgcolor="rgba(255,255,255,255)",
        #                     plot_bgcolor="rgba(255, 255, 255, 255)",
        #                     yaxis=dict(showline=True, linewidth=1, linecolor='black'),
        #                     xaxis=dict(showline=True, linewidth=1, linecolor='black')))
        #     fig.update_xaxes(categoryorder='array', categoryarray = plot_options["x_axis_order"])

        # el
        if plot_options["ID mode"] == "stacked":
            
            
            last_layer = None
            layers = plot_options["cutoffs"]
            fig_data = []
            i = 0
            for eachLayer in layers: 
                if last_layer == None:
                    fig_data.append(go.Bar(
                        name = eachLayer,
                        x = ID_data.loc[(ID_data["Cutoff"]==eachLayer),plot_options["Group By X"]].tolist(),
                        y = ID_data.loc[(ID_data["Cutoff"]==eachLayer),"IDs"].tolist(),
                        marker_color = plot_options["color"][i],
                        text = [round(x) for x in ID_data.loc[(ID_data["Cutoff"]==eachLayer),"IDs"].tolist()],
                        
                    ))
                    bases = ID_data.loc[(ID_data["Cutoff"]==eachLayer),"IDs"]
                    i = i + 1
                    
                else:
                    fig_data.append(go.Bar(
                        name = eachLayer,
                        x = ID_data.loc[(ID_data["Cutoff"]==eachLayer),plot_options["Group By X"]].tolist(),
                        y = ID_data.loc[(ID_data["Cutoff"]==eachLayer),"IDs"].tolist(),
                        base=bases,
                        marker_color = plot_options["color"][i],
                        text = [round(x) for x in bases + ID_data.loc[(ID_data["Cutoff"]==eachLayer),"IDs"].tolist()],
                    ))
                    # print(bases)
                    bases = bases + ID_data.loc[(ID_data["Cutoff"]==eachLayer),"IDs"].tolist()
                i = i + 1
                last_layer = eachLayer
            fig = go.Figure(
                    data = fig_data,
                    layout=go.Layout(
                    yaxis_title=plot_options["Y Title"],
                    xaxis_title=plot_options["Group By X"],
                    barmode="stack", 
                    paper_bgcolor="rgba(255,255,255,255)",
                    plot_bgcolor="rgba(255, 255, 255, 255)",
                    yaxis=dict(showline=True, linewidth=1, linecolor='black'),
                    xaxis=dict(showline=True, linewidth=1, linecolor='black')
                ))          
            fig.update_xaxes(categoryorder='array', categoryarray = plot_options["x_axis_order"])            
        # elif plot_options["ID mode"] == "grouped_stacked":
        #     # plot options
        #     # error bar
        #     if plot_options["error bar"] == "stdev":
        #         error_bars = "stdev"
        #         error_visibile = True
        #     elif plot_options["error bar"] == "ci95":
        #         error_bars = "confInt"
        #         error_visibile = True
        #     else:
        #         error_bars = "stdev"
        #         error_visibile = False

        #     # mean label
        #     if plot_options["mean label"] == "True" or \
        #             plot_options["mean label"] == True:
        #         total_labels = [{"x": x, "y": total*1.15, "text": str(
        #             int(total)), "showarrow": False} for x, total in zip(
        #                 ID_data["Conditions"], ID_data["IDs"])]
        #     else:
        #         total_labels = []   # no mean labels

        #     if plot_options["Group By X"] == "ID_Mode" or plot_options["Group By Color"] == "ID_Mode"or plot_options["Group By Stack"] == "ID_Mode":  # total separated
        #         ID_data = ID_data[ID_data["ID_Mode"] != "Total_IDs"]
        #     else:
        #         ID_data = ID_data[ID_data["ID_Mode"] == "Total_IDs"]

        #     #make data tidy
        #     if plot_options["Group By Stack"] == "ID_Mode":
        #         layers = ["MS2_IDs", "MBR_IDs"]
        #     else:
        #         layers = ID_data.groupby(plot_options["Group By Stack"]).first().reset_index()[plot_options["Group By Stack"]].tolist()
        #     categories = plot_options["color_order"]
            
        #     fig_data = []
        #     i = 0
        #     for eachCategory in categories:
        #         last_layer = None
        #         j = 0
        #         for eachLayer in layers: 
        #             if last_layer == None:
        #                 fig_data.append(go.Bar(
        #                     name = str(eachLayer) + " " + str(eachCategory),
        #                     x = ID_data.loc[(ID_data[plot_options["Group By Color"]]==eachCategory)&(ID_data[plot_options["Group By Stack"]]==eachLayer),plot_options["Group By X"]],
        #                     y = ID_data.loc[(ID_data[plot_options["Group By Color"]]==eachCategory)&(ID_data[plot_options["Group By Stack"]]==eachLayer),"IDs"],
        #                     offsetgroup=i,
        #                     text = [round(x) for x in ID_data.loc[(ID_data[plot_options["Group By Color"]]==eachCategory)&(ID_data[plot_options["Group By Stack"]]==eachLayer),"IDs"].tolist()],
        #                     marker_color = plot_options["color"][i],
        #                     error_y = dict(
        #                         type = "data",
        #                         array = ID_data.loc[(ID_data[plot_options["Group By Color"]]==eachCategory)&(ID_data[plot_options["Group By Stack"]]==eachLayer),error_bars],
        #                         visible=True)
        #                 ))
        #                 bases = ID_data.loc[(ID_data[plot_options["Group By Color"]]==eachCategory)&(ID_data[plot_options["Group By Stack"]]==eachLayer),"IDs"]
        #             else:

        #                 fig_data.append(go.Bar(
        #                     name = str(eachLayer) + " " + str(eachCategory),
        #                     x = ID_data.loc[(ID_data[plot_options["Group By Color"]]==eachCategory)&(ID_data[plot_options["Group By Stack"]]==eachLayer),plot_options["Group By X"]],
        #                     y = ID_data.loc[(ID_data[plot_options["Group By Color"]]==eachCategory)&(ID_data[plot_options["Group By Stack"]]==eachLayer),"IDs"],
        #                     base=bases,
        #                     offsetgroup=i,
        #                     text = [round(x) for x in ID_data.loc[(ID_data[plot_options["Group By Color"]]==eachCategory)&(ID_data[plot_options["Group By Stack"]]==eachLayer),"IDs"].tolist()+bases],
        #                     marker_color = plot_options["color"][i],
        #                     opacity=1/2**j,
        #                     error_y = dict(
        #                         type = "data",
        #                         array = ID_data.loc[(ID_data[plot_options["Group By Color"]]==eachCategory)&(ID_data[plot_options["Group By Stack"]]==eachLayer),error_bars],
        #                         visible=True)
        #                     ))
        #                 bases = bases + ID_data.loc[(ID_data[plot_options["Group By Color"]]==eachCategory)&(ID_data[plot_options["Group By Stack"]]==eachLayer),"IDs"].tolist()
        #             last_layer=eachLayer
        #             j = j + 1
        #         i = i + 1

        #     fig = go.Figure(
        #         fig_data,
        #         layout=go.Layout(
        #             yaxis_title=plot_options["Y Title"],
        #             xaxis_title=plot_options["Group By X"],
        #             barmode="group",
        #             plot_bgcolor="rgba(255, 255, 255, 255)",
        #             paper_bgcolor="rgba(255, 255, 255, 255)",
        #             yaxis=dict(showline=True, linewidth=1, linecolor='black'),
        #             xaxis=dict(showline=True, linewidth=1, linecolor='black')
        #         )
        #     )
        #     fig.update_xaxes(categoryorder='array', categoryarray = plot_options["x_axis_order"])
        # else:

        #     # plot options
        #     # error bar
        #     if plot_options["error bar"] == "stdev":
        #         error_bars = "stdev"
        #     elif plot_options["error bar"] == "ci95":
        #         error_bars = "confInt"
        #     else:
        #         error_bars = None

        #     # mean label
        #     if plot_options["mean label"] == "True" or \
        #             plot_options["mean label"] == True:
        #         total_labels = [{"x": x, "y": total*1.15, "text": str(
        #             int(total)), "showarrow": False} for x, total in zip(
        #                 ID_data["Conditions"], ID_data["IDs"])]
        #     else:
        #         total_labels = []   # no mean labels

        #     # create the plot
        #     fig = px.bar(ID_data,
        #                  x="Conditions",
        #                  y="IDs",
        #                  error_y=error_bars,
        #                  color="Conditions",
        #                  color_discrete_sequence=plot_options["color"],
        #                  width=plot_options["width"],
        #                  height=plot_options["height"],
        #                  )
        #     fig.update_layout(xaxis_title=plot_options["X Title"],
        #                       yaxis_title=plot_options["Y Title"],
        #                       annotations=total_labels,
        #                       font=plot_options["font"],
        #                       plot_bgcolor ="rgba(255, 255, 255, 255)",
        #                       paper_bgcolor="rgba(255, 255, 255, 255)",
        #                       yaxis=dict(showline=True, linewidth=1, linecolor='black'),
        #                       xaxis=dict(showline=True, linewidth=1, linecolor='black')
        #                       )
        #     fig.update_xaxes(categoryorder='array', categoryarray = plot_options["x_axis_order"])
        return fig

    # CV Violin plots ###
    def CV_plots(self, data_object, plot_options, saved_settings, username=None):
        """_Prepare data for creating protein CV violin plots_
        """
        group_names = [key for key in saved_settings.keys() if "Order@" not in str(key)]
        
        # import the data and save order
        group_dict = {}

        if plot_options["CV mode"] == "MS2" or plot_options["CV mode"] == "total" or plot_options["CV mode"] == "stacked":
            x_axis_order = saved_settings["Order@Conditions"]
        else:
            x_axis_order = saved_settings["Order@"+plot_options["Group By X"]]
        plot_options["x_axis_order"] = x_axis_order
        if plot_options["CV mode"] == "grouped" or plot_options["CV mode"] == "grouped_stacked" and plot_options["Group By Color"] != "ID_Mode":
            color_order = saved_settings["Order@"+plot_options["Group By Color"]]
            plot_options["color_order"] = color_order
        # filter runs into different groups
        i = 1
        runname_list = []  # contain list of run names list for each groups
        for eachGroup in group_names:

            run_ids = saved_settings[eachGroup]["records"]
            runname_sublist = data_object["run_metadata"].loc[data_object["run_metadata"]["Run Identifier"].isin(run_ids)]["Run Names"].to_list()
            if self.processor.data_type == "TMT":
                channel_ids = data_object["run_metadata"].loc[data_object["run_metadata"]["Run Identifier"].isin(run_ids)]["Channel Identifier"].to_list()
                group_dict[eachGroup] =  self.processor.filter_by_channel_id(
                    data_object,
                    channel_ids)  
            elif self.processor.data_type == "LF":

                group_dict[eachGroup] =  self.processor.filter_by_id(
                    data_object,
                    run_ids)  # prevent the list from being changed
            runname_list.append(runname_sublist)
            i += 1

        if plot_options["plot_type"] == 1:
            matrix_name = "protein_abundance"
        elif plot_options["plot_type"] == 2:
            matrix_name = "peptide_abundance" 

        # create a dictionary to store the intensity data
        Intensity_dict = {}
        if plot_options["Group By X"] == "ID_Mode"or plot_options["Group By Color"] == "ID_Mode" \
            or plot_options["Group By Stack"] == "ID_Mode"and not (plot_options["CV mode"] == "total" or plot_options["CV mode"] == "MS2"):  # total separated
            Intensity_dict_MS2 = {}
            for eachGroup in group_names:
                current_condition_data =  self.processor.filter_by_missing_values(
                    group_dict[eachGroup])
                Intensity_dict[eachGroup] =  self.processor.NormalizeToMedian(
                    current_condition_data[matrix_name])
                current_condition_data_MS2 =  self.processor.filter_by_missing_values_MS2( #returns
                    group_dict[eachGroup])
                Intensity_dict_MS2[eachGroup] =  self.processor.NormalizeToMedian(
                    current_condition_data_MS2[matrix_name])
        elif plot_options["CV mode"] == "MS2":
            Intensity_dict_MS2 = {}
            for eachGroup in group_names:
                current_condition_data_MS2 =  self.processor.filter_by_missing_values_MS2(
                    group_dict[eachGroup])
                Intensity_dict[eachGroup] =  self.processor.NormalizeToMedian(
                    current_condition_data_MS2[matrix_name])
        else:
            for eachGroup in group_names:
                current_condition_data =  self.processor.filter_by_missing_values(
                    group_dict[eachGroup])
                Intensity_dict[eachGroup] =  self.processor.NormalizeToMedian(
                    current_condition_data[matrix_name])

        all_cvs = pd.DataFrame()

        for eachGroup in Intensity_dict:
            current =  self.processor.calculate_cvs(
                Intensity_dict[eachGroup]).assign(Conditions=eachGroup,ID_Mode="All IDs")
            all_cvs = pd.concat([all_cvs, current], ignore_index=True)
        if plot_options["Group By X"] == "ID_Mode"or plot_options["Group By Color"] == "ID_Mode"or plot_options["Group By Stack"] == "ID_Mode"and not plot_options["CV mode"] == "total":  # total separated  # total separated
            print("MS2 CVs...")
            for eachGroup in Intensity_dict_MS2:
                current =  self.processor.calculate_cvs(
                    Intensity_dict_MS2[eachGroup]).assign(Conditions=eachGroup,ID_Mode="MS2 IDs")
                    
                all_cvs = pd.concat([all_cvs, current], ignore_index=True)

        #add columns for the categories specified in settings file (the one with all the filenames)
        standard_groups = ["filter_in","filter_out","records"]
        categories = [col for col in list(saved_settings[list(group_names)[0]].keys()) if col not in standard_groups]
        for eachCategory in categories:
            all_cvs[eachCategory] = ""
            for eachGroup in group_names:
                all_cvs.loc[all_cvs["Conditions"]==eachGroup,eachCategory] = saved_settings[eachGroup][eachCategory]


        # ######################all_CVs format###################
    #      Accession     intensity          stdev          CV   Conditions
    # 0       A6NHR9  3.248547e+06  672989.819300   20.716643    DDMandDTT
    # 1       A8MTJ3  5.031539e+05  195535.383583   38.861944    DDMandDTT
    # 2       E9PAV3  5.330290e+05  161385.491163   30.277056    DDMandDT
        #######################################################

        return self.plot_CV_violin(allCVs=all_cvs,
                            username=username,
                            plot_options=plot_options,
                            saved_settings=saved_settings)


    def plot_CV_violin(self, allCVs, 
                    username=None,
                    plot_options=None,
                    saved_settings=None
                    ):
        """_Plot the CV violin plot for the given data._

        Args:
            allCVs (_type_): _description_
            username (_type_, optional): _description_. Defaults to None.
            plot_options (_type_, optional): _description_. Defaults to None.

        Returns:
            _type_: _description_
        """
        CSV_link = None
        SVG_link = None

        group_names = [key for key in saved_settings.keys() if "Order@" not in str(key)]

        allCVs_summary = allCVs.groupby(["Conditions"]).agg(
            {'CV': ['median', 'mean']}).reset_index()
        allCVs_summary["ID_Mode"] = "All IDs"
        temp = allCVs[allCVs["ID_Mode"]=="MS2 IDs"].groupby(["Conditions"]).agg(
            {'CV': ['median', 'mean']}).reset_index()
        temp["ID_Mode"] = "MS2 IDs"
        allCVs_summary = pd.concat([temp,allCVs_summary])
        allCVs_summary.columns = ["Conditions", 'meds', 'CoVar',"ID_Mode"]

        
        standard_groups = ["filter_in","filter_out","records"]
        categories = [col for col in list(saved_settings[group_names[0]].keys()) if col not in standard_groups]
        for eachCategory in categories:
            allCVs_summary[eachCategory] = ""
            for eachGroup in group_names:
                allCVs_summary.loc[allCVs_summary["Conditions"]==eachGroup,eachCategory] = saved_settings[eachGroup][eachCategory]


        if plot_options["CV mode"] == "grouped":  
            #find out present categories
            categories = plot_options["color_order"] 
            # create the plot
            fig_data = []
            #display(ID_data)
            i = 0
            for eachCategory in categories:
                fig_data.append(go.Violin(name = str(eachCategory),
                            x=allCVs.loc[allCVs[plot_options["Group By Color"]]==eachCategory,
                            plot_options["Group By X"]].tolist(),
                            y=allCVs.loc[allCVs[plot_options["Group By Color"]]==eachCategory,"CV"].tolist(),
                            fillcolor = plot_options["color"][i],                 
                            box=dict(visible=bool(plot_options["box"]))
                            ))
                i = i + 1
            
                            
            fig = go.Figure(data = fig_data)
            
            fig.update_layout(violinmode='group')
            fig.update_xaxes(categoryorder='array', categoryarray = plot_options["x_axis_order"])
        elif plot_options["CV mode"] == "stacked":  
                
            #find out present categories
            layers = allCVs.groupby(plot_options["Group By Stack"]).first().reset_index()[plot_options["Group By Stack"]].tolist()
            # create the plot
            fig_data = []
            #display(ID_data)
            i = 0

            for eachLayer in layers:
                fig_data.append(go.Violin(name = eachLayer,
                            x=allCVs.loc[allCVs[plot_options["Group By Stack"]]==eachLayer,plot_options["Group By X"]].tolist(),
                            y=allCVs.loc[allCVs[plot_options["Group By Stack"]]==eachLayer,"CV"].tolist(),
                            box=dict(visible=bool(plot_options["box"])),
                            fillcolor = plot_options["color"][i]
                            ))
                i = i + 1
            fig = go.Figure(data = fig_data)
            
            fig.update_layout(violinmode='overlay',
                            plot_bgcolor='white',
                            paper_bgcolor='white',
                            yaxis=dict(showline=True, linewidth=1, linecolor='black'),
                            xaxis=dict(showline=True, linewidth=1, linecolor='black')
            )  
            fig.update_xaxes(categoryorder='array', categoryarray = plot_options["x_axis_order"])

        elif plot_options["CV mode"] == "grouped_stacked":
            #make data tidy
            layers = allCVs.groupby(plot_options["Group By Stack"]).first().reset_index()[plot_options["Group By Stack"]].tolist()
            categories = plot_options["color_order"] 
            fig_data = []
            i = 0
            j = 0
            for eachCategory in categories:
                for eachLayer in layers:
                    fig_data.append(go.Violin(
                        name = str(eachLayer) + " " + str(eachCategory),
                        x = allCVs.loc[(allCVs[plot_options["Group By Color"]]==eachCategory)&(allCVs[plot_options["Group By Stack"]]==eachLayer),plot_options["Group By X"]],
                        y = allCVs.loc[(allCVs[plot_options["Group By Color"]]==eachCategory)&(allCVs[plot_options["Group By Stack"]]==eachLayer),"CV"],
                        box=dict(visible=bool(plot_options["box"])),
                        offsetgroup=i,
                        fillcolor = plot_options["color"][j]
                        ))
                    j = j + 1
                i = i + 1

            fig = go.Figure(
                fig_data,
                layout=go.Layout(
                    yaxis_title=plot_options["Y Title"],
                    xaxis_title=plot_options["Group By X"],
                    violinmode="group",
                    plot_bgcolor='white',
                    paper_bgcolor='white',
                    yaxis=dict(showline=True, linewidth=1, linecolor='black'),
                    xaxis=dict(showline=True, linewidth=1, linecolor='black')
                )
            )
            fig.update_xaxes(categoryorder='array', categoryarray = plot_options["x_axis_order"])
        else:
        # create the interactive plot
            # median label
            if plot_options["median label"] == "True" or \
                    plot_options["median label"] == True:
                total_labels = [{"x": x, "y": total*1.15, "text": str(
                    round(total,1)), "showarrow": False} for x, total in zip(
                    allCVs_summary["Conditions"], allCVs_summary["meds"])]
            else:
                total_labels = []   # no median labels
            fig = px.violin(allCVs,
                            x="Conditions",
                            y='CV',
                            color="Conditions",
                            box=bool(plot_options["box"]),
                            hover_data=["Conditions", 'CV'],
                            color_discrete_sequence=plot_options["color"],
                            width=plot_options["width"],
                            height=plot_options["height"],
                            )

            fig.update_layout(
                yaxis=dict(title=plot_options["Y Title"],
                        range=plot_options["ylimits"], showline=True, linewidth=1, linecolor='black'),
                font=plot_options["font"],
                xaxis=dict(title=plot_options["X Title"], showline=True, linewidth=1, linecolor='black'),
                showlegend=True,
                annotations=total_labels,
                plot_bgcolor='white',
                paper_bgcolor='white',
                )
            fig.update_xaxes(categoryorder='array', categoryarray = plot_options["x_axis_order"])

        if self.write_output:        
            # create the file for donwnload
            img_dir = os.path.join(self.app_folder, "images/")
            if not os.path.exists(img_dir):
                Path(img_dir).mkdir(parents=True)

            fig.write_image(os.path.join(
                img_dir, f"{username}_CV_Violin_Plot.svg"), format = "svg", validate = False, engine = "kaleido")
            
            # create the download CSV and its link
            data_dir = os.path.join(self.app_folder, "csv/")
            if not os.path.exists(data_dir):
                Path(data_dir).mkdir(parents=True)
            allCVs.to_csv(os.path.join(
                data_dir, f"{username}_all_CV.csv"), index=False)
            allCVs_summary.to_csv(os.path.join(
                data_dir, f"{username}_CV_summary.csv"), index=False)
            print("Downloading links...")
            CSV_link = f"/files/{self.url_base}/csv/" \
                f"{username}_all_CV.csv"

            # download SVG link
            SVG_link = f"/files/{self.url_base}/images/" \
                f"{username}_CV_Violin_Plot.svg"


        return fig, CSV_link, SVG_link


    def inclusion_venn_plots(self, data_object, plot_options, saved_settings, username=None,miss_val_thresh=33):
        """_Prepare data for creating ID veens plots (up to three groups)_
        """
        group_names = []

        # no compare groups is provided, compare first two
        
        for each_key in plot_options["compare groups"]:
            if each_key in plot_options["compare groups"] and saved_settings["Order@Conditions"]:
                group_names.append(each_key)
        # import the data
        group_dict = {}

        # filter runs into different groups
        i = 0
        runname_list = []  # contain list of run names list for each groups
        for eachGroup in group_names:

            run_ids = saved_settings[eachGroup]["records"]
            runname_sublist = data_object["run_metadata"].loc[data_object["run_metadata"]["Run Identifier"].isin(run_ids)]["Run Names"].to_list()
            
            if self.processor.data_type == "TMT":
                channel_ids = data_object["run_metadata"].loc[data_object["run_metadata"]["Run Identifier"].isin(run_ids)]["Channel Identifier"].to_list()
                group_dict[eachGroup] =  self.processor.filter_by_channel_id(
                    data_object,
                    channel_ids)  
            elif self.processor.data_type == "LF":

                group_dict[eachGroup] =  self.processor.filter_by_id(
                    data_object,
                    run_ids) # prevent the list from being changed
            runname_list.append(runname_sublist)
            i += 1
        
        labels_set = ["Inclusion List"]
        Inclusion_List = pd.read_table(plot_options["inclusion list"],sep=",")["Accession"].to_list()
        data_set = [set(Inclusion_List)]
        if plot_options["plot_type"] == 1:
            matrix_name = "protein_abundance"
            molecule_name = "Accession"
            is_protein = True
        elif plot_options["plot_type"] == 2:
            matrix_name = "peptide_abundance"
            molecule_name = "Annotated Sequence"
            is_protein = False
        
        for eachGroup in group_names:
            
            current_condition_data = self.processor.filter_by_missing_values_MS2(
                group_dict[eachGroup], is_protein=is_protein, missing_value_thresh=miss_val_thresh)

            data_set.append(
                set(current_condition_data[matrix_name][molecule_name].unique()))
            labels_set.append(eachGroup)

        #print(data_set)

        fig = self.venn_to_plotly(
            data_set,
            labels_set,
            plot_options=plot_options,
            username=username)
        CSV_link = None
        SVG_link = None

        if self.write_output:
            print("Downloading links...")
            # SVG file link
            SVG_link = f"/files/{self.url_base}/images/" \
                f"{username}_ID_venns_Plot.svg"

            # create the file for donwnload
            img_dir = os.path.join(self.app_folder, "images/")
            if not os.path.exists(img_dir):
                Path(img_dir).mkdir(parents=True)

            fig.write_image(os.path.join(
                img_dir, f"{username}_ID_venns_Plot.svg"), format = "svg", validate = False, engine = "kaleido")
            
            data_dir = os.path.join(self.app_folder, "csv/")
            if not os.path.exists(data_dir):
                Path(data_dir).mkdir(parents=True)
            i = 0
            for eachSet in data_set:
                pd.DataFrame({"Accession": list(eachSet)}).to_csv(os.path.join(
                    data_dir, f"{username}_{labels_set[i]}_Venn.csv"), index=False)
                i = i + 1
        
        return fig, SVG_link, CSV_link

    def venn_to_plotly(self, L_sets,
                    L_labels=None,
                    plot_options=None,
                    username=None):
        """_Creates a venn diagramm from a list of
        sets and returns a plotly figure_
        """
        
        # get number of sets
        n_sets = len(L_sets)

        # choose and create matplotlib venn diagramm
        if n_sets == 2:
            if L_labels and len(L_labels) == n_sets:
                v = venn2(L_sets, L_labels)
            else:
                v = venn2(L_sets)
        elif n_sets == 3:
            if L_labels and len(L_labels) == n_sets:
                v = venn3(L_sets, L_labels)
            else:
                v = venn3(L_sets)
        # supress output of venn diagramm
        # plt.show()
        plt.close()

        # Create empty lists to hold shapes and annotations
        L_shapes = []
        L_annotation = []

        # Define color list for sets
        L_color = plot_options["color"]

        # Create empty list to make hold of min and max values of set shapes
        L_x_max = []
        L_y_max = []
        L_x_min = []
        L_y_min = []

        for i in range(0, n_sets):

            # create circle shape for current set

            shape = go.layout.Shape(
                type="circle",
                xref="x",
                yref="y",
                x0=v.centers[i][0] - v.radii[i],
                y0=v.centers[i][1] - v.radii[i],
                x1=v.centers[i][0] + v.radii[i],
                y1=v.centers[i][1] + v.radii[i],
                fillcolor=L_color[i],
                line_color=L_color[i],
                opacity=plot_options["opacity"]
            )

            L_shapes.append(shape)

            # create set label for current set
            try:
                anno_set_label = go.layout.Annotation(
                    xref="x",
                    yref="y",
                    x=v.set_labels[i].get_position()[0],
                    y=v.set_labels[i].get_position()[1],
                    text=v.set_labels[i].get_text(),
                    showarrow=False
                )

                L_annotation.append(anno_set_label)

                # get min and max values of current set shape
                L_x_max.append(v.centers[i][0] + v.radii[i])
                L_x_min.append(v.centers[i][0] - v.radii[i])
                L_y_max.append(v.centers[i][1] + v.radii[i])
                L_y_min.append(v.centers[i][1] - v.radii[i])
            except Exception as err:
                print(f"No set labels found {err}")

        # determine number of subsets
        n_subsets = sum([scipy.special.binom(n_sets, i+1)
                        for i in range(0, n_sets)])

        for i in range(0, int(n_subsets)):
            try:

                # create subset label (number of common elements for current subset

                anno_subset_label = go.layout.Annotation(
                    xref="x",
                    yref="y",
                    x=v.subset_labels[i].get_position()[0],
                    y=v.subset_labels[i].get_position()[1],
                    text=v.subset_labels[i].get_text(),
                    showarrow=False
                )

                L_annotation.append(anno_subset_label)
            except Exception as err:
                print(f"No set labels found {err}")
        # define off_set for the figure range
        off_set = 0.2

        # get min and max for x and y dimension to set the figure range
        x_max = max(L_x_max) + off_set
        x_min = min(L_x_min) - off_set
        y_max = max(L_y_max) + off_set
        y_min = min(L_y_min) - off_set

        # create plotly figure

        fig = go.Figure()

        # set xaxes range and hide ticks and ticklabels
        fig.update_xaxes(
            range=[x_min, x_max],
            showticklabels=False,
            ticklen=0
        )

        # set yaxes range and hide ticks and ticklabels
        fig.update_yaxes(
            range=[y_min, y_max],
            scaleanchor="x",
            scaleratio=1,
            showticklabels=False,
            ticklen=0
        )

        # set figure properties and add shapes and annotations
        fig.update_layout(
            plot_bgcolor='white',
            margin=dict(b=0, l=10, pad=0, r=10, t=40),
            width=800,
            height=400,
            shapes=L_shapes,
            annotations=L_annotation,
            title=dict(text=plot_options["title"], x=0.5, xanchor='center')
        )
        

        return fig
    
    def venns_plots(self, data_object, plot_options, saved_settings, username=None,miss_val_thresh=33):
        """_Prepare data for creating ID veens plots (up to three groups)_
        """
        group_names = []

        # no compare groups is provided, compare first two
        
        for each_key in plot_options["compare groups"]:
            if each_key in plot_options["compare groups"] and saved_settings["Order@Conditions"]:
                group_names.append(each_key)
        # import the data
        group_dict = {}

        # filter runs into different groups
        i = 0
        runname_list = []  # contain list of run names list for each groups
        for eachGroup in group_names:

            run_ids = saved_settings[eachGroup]["records"]
            runname_sublist = data_object["run_metadata"].loc[data_object["run_metadata"]["Run Identifier"].isin(run_ids)]["Run Names"].to_list()
            
            if self.processor.data_type == "TMT":
                print(data_object["protein_ID_matrix"])
                channel_ids = data_object["run_metadata"].loc[data_object["run_metadata"]["Run Identifier"].isin(run_ids)]["Channel Identifier"].to_list()
                group_dict[eachGroup] =  self.processor.filter_by_channel_id(
                    data_object,
                    channel_ids)  
                print(group_dict[eachGroup]["protein_ID_matrix"])
            elif self.processor.data_type == "LF":

                group_dict[eachGroup] =  self.processor.filter_by_id(
                    data_object,
                    run_ids) # prevent the list from being changed
            runname_list.append(runname_sublist)
            i += 1
        data_set = []
        labels_set = []
        if plot_options["plot_type"] == 1:
            matrix_name = "protein_abundance"
            molecule_name = "Accession"
            is_protein = True
        elif plot_options["plot_type"] == 2:
            matrix_name = "peptide_abundance"
            molecule_name = "Annotated Sequence"
            is_protein = False
        
        for eachGroup in group_names:
            
            current_condition_data = self.processor.filter_by_missing_values_MS2(
                group_dict[eachGroup], is_protein=is_protein, missing_value_thresh=miss_val_thresh)

            data_set.append(
                set(current_condition_data[matrix_name][molecule_name].unique()))
            labels_set.append(eachGroup)

        #print(data_set)

        fig = self.venn_to_plotly(
            data_set,
            labels_set,
            plot_options=plot_options,
            username=username)
        CSV_link = None
        SVG_link = None

        if self.write_output:
            print("Downloading links...")
            # SVG file link
            SVG_link = f"/files/{self.url_base}/images/" \
                f"{username}_ID_venns_Plot.svg"

            # create the file for donwnload
            img_dir = os.path.join(self.app_folder, "images/")
            if not os.path.exists(img_dir):
                Path(img_dir).mkdir(parents=True)

            fig.write_image(os.path.join(
                img_dir, f"{username}_ID_venns_Plot.svg"), format = "svg", validate = False, engine = "kaleido")
            
            data_dir = os.path.join(self.app_folder, "csv/")
            if not os.path.exists(data_dir):
                Path(data_dir).mkdir(parents=True)
            i = 0
            for eachSet in data_set:
                pd.DataFrame({"Accession": list(eachSet)}).to_csv(os.path.join(
                    data_dir, f"{username}_{labels_set[i]}_Venn.csv"), index=False)
                i = i + 1
        
        return fig, SVG_link, CSV_link

    def venn_to_plotly(self, L_sets,
                    L_labels=None,
                    plot_options=None,
                    username=None):
        """_Creates a venn diagramm from a list of
        sets and returns a plotly figure_
        """
        
        # get number of sets
        n_sets = len(L_sets)

        # choose and create matplotlib venn diagramm
        if n_sets == 2:
            if L_labels and len(L_labels) == n_sets:
                v = venn2(L_sets, L_labels)
            else:
                v = venn2(L_sets)
        elif n_sets == 3:
            if L_labels and len(L_labels) == n_sets:
                v = venn3(L_sets, L_labels)
            else:
                v = venn3(L_sets)
        # supress output of venn diagramm
        # plt.show()
        plt.close()

        # Create empty lists to hold shapes and annotations
        L_shapes = []
        L_annotation = []

        # Define color list for sets
        L_color = plot_options["color"]

        # Create empty list to make hold of min and max values of set shapes
        L_x_max = []
        L_y_max = []
        L_x_min = []
        L_y_min = []

        for i in range(0, n_sets):

            # create circle shape for current set

            shape = go.layout.Shape(
                type="circle",
                xref="x",
                yref="y",
                x0=v.centers[i][0] - v.radii[i],
                y0=v.centers[i][1] - v.radii[i],
                x1=v.centers[i][0] + v.radii[i],
                y1=v.centers[i][1] + v.radii[i],
                fillcolor=L_color[i],
                line_color=L_color[i],
                opacity=plot_options["opacity"]
            )

            L_shapes.append(shape)

            # create set label for current set
            try:
                anno_set_label = go.layout.Annotation(
                    xref="x",
                    yref="y",
                    x=v.set_labels[i].get_position()[0],
                    y=v.set_labels[i].get_position()[1],
                    text=v.set_labels[i].get_text(),
                    showarrow=False
                )

                L_annotation.append(anno_set_label)

                # get min and max values of current set shape
                L_x_max.append(v.centers[i][0] + v.radii[i])
                L_x_min.append(v.centers[i][0] - v.radii[i])
                L_y_max.append(v.centers[i][1] + v.radii[i])
                L_y_min.append(v.centers[i][1] - v.radii[i])
            except Exception as err:
                print(f"No set labels found {err}")

        # determine number of subsets
        n_subsets = sum([scipy.special.binom(n_sets, i+1)
                        for i in range(0, n_sets)])

        for i in range(0, int(n_subsets)):
            try:

                # create subset label (number of common elements for current subset

                anno_subset_label = go.layout.Annotation(
                    xref="x",
                    yref="y",
                    x=v.subset_labels[i].get_position()[0],
                    y=v.subset_labels[i].get_position()[1],
                    text=v.subset_labels[i].get_text(),
                    showarrow=False
                )

                L_annotation.append(anno_subset_label)
            except Exception as err:
                print(f"No set labels found {err}")
        # define off_set for the figure range
        off_set = 0.2

        # get min and max for x and y dimension to set the figure range
        x_max = max(L_x_max) + off_set
        x_min = min(L_x_min) - off_set
        y_max = max(L_y_max) + off_set
        y_min = min(L_y_min) - off_set

        # create plotly figure

        fig = go.Figure()

        # set xaxes range and hide ticks and ticklabels
        fig.update_xaxes(
            range=[x_min, x_max],
            showticklabels=False,
            ticklen=0
        )

        # set yaxes range and hide ticks and ticklabels
        fig.update_yaxes(
            range=[y_min, y_max],
            scaleanchor="x",
            scaleratio=1,
            showticklabels=False,
            ticklen=0
        )

        # set figure properties and add shapes and annotations
        fig.update_layout(
            plot_bgcolor='white',
            margin=dict(b=0, l=10, pad=0, r=10, t=40),
            width=800,
            height=400,
            shapes=L_shapes,
            annotations=L_annotation,
            title=dict(text=plot_options["title"], x=0.5, xanchor='center')
        )
        

        return fig
    # ###HYE plots####
    def HYE_plots(self, data_object,  plot_options, saved_settings, username=None):
        """_Prepare data for creating intensity volcano plots (two groups)_
        """
        group_names = []

        # no compare groups is provided, compare first two
        for each_key in saved_settings["Order@Conditions"]:
                group_names.append(each_key)

        # import the data
        group_dict = {}

        # filter runs into different groups
        i = 0
        runname_list = []  # contain list of run names list for each groups
        for eachGroup in group_names:
            runname_sublist = saved_settings[eachGroup]["records"]

            if self.processor.data_type == "TMT":
                channel_ids = data_object["run_metadata"].loc[data_object["run_metadata"]["Run Identifier"].isin(runname_sublist)]["Channel Identifier"].to_list()
                group_dict[eachGroup] =  self.processor.filter_by_channel_id(
                    data_object,
                    channel_ids)  
            elif self.processor.data_type == "LF":

                group_dict[eachGroup] =  self.processor.filter_by_id(
                    data_object,
                    runname_sublist)  # prevent the list from being changed
            runname_list.append(runname_sublist)
            i += 1
        # create a dictionary to store the intensity data
        Intensity_dict = {}

        for eachGroup in group_names:
            current_condition_data = self.processor.filter_by_missing_values(
                group_dict[eachGroup])
            Intensity_dict[eachGroup] = self.processor.NormalizeToMedian(
                current_condition_data["protein_abundance"],apply_log2=True)
            
                
            
            # calculate mean, standard deviation, and the number of non-null
            # elements for each row/protein
            Intensity_dict[eachGroup][eachGroup+ "_Intensity"] = Intensity_dict[eachGroup].iloc[1:].mean(axis=1)
            Intensity_dict[eachGroup][eachGroup+ "_stdev"]= Intensity_dict[eachGroup].iloc[1:-1].mean(axis=1)
            Intensity_dict[eachGroup][eachGroup+ "_count"]= Intensity_dict[eachGroup].iloc[1:-2].shape[1]-Intensity_dict[eachGroup].iloc[1:-2].isna().sum(axis=1)
            Intensity_dict[eachGroup] = Intensity_dict[eachGroup].loc[:,["Accession",
                                                                        eachGroup+"_Intensity",
                                                                        eachGroup+"_stdev",
                                                                        eachGroup+"_count",
                                                                        ]]
            # print(Intensity_dict[eachGroup])
                
            """ group1Data
                    group1_Intensity  group1_stdev  group1_num   Accession
            0        2.824766e+05  1.708060e+05          15  A0A0B4J2D5
            1        2.650998e+06  6.259645e+05          15      A2RUR9
            2        1.973150e+05  5.645698e+04          15      A8MTJ3
            3        2.524020e+05  1.355699e+05          15      A8MWD9
            """

        #this is the median for each protein across both groups, you could do a ratio, but this puts it in terms of
        # allmedian = pd.DataFrame(Intensity_dict[group][group+'_Intensity'] for group in list(Intensity_dict.keys())).median(axis=1,numeric_only=True)


        # commonProts=pd.DataFrame()
        # for eachGroup in group_names:
        #     these_prots = Intensity_dict[eachGroup].loc[:, ['Accession']]

        #     if commonProts.empty:
        #         commonProts = these_prots
        #     # find common proteins
        #     else:
        #         commonProts = (commonProts.merge(these_prots, on='Accession', how='inner'))
        #     # only leave common proteins
            

        # volcanoData = pd.DataFrame()

        # for eachGroup in group_names:
        #     # Intensity_dict[eachGroup] = (Intensity_dict[eachGroup].merge(commonProts, on='Accession', how='inner'))
        #     # print(Intensity_dict[eachGroup])
        #     # this_median = Intensity_dict[eachGroup][eachGroup+'_Intensity'].median(axis=1,numeric_only =True)  
        #     # # calculate the ratio between two group median,
        #     # # will be used to normalize them
        #     # this_ratio = allmedian / this_median

        #     # merge these two set of data together, adjust groups with ratio to 
        #     # median of all. 
            
        #     # Intensity_dict[eachGroup][eachGroup+"_Intensity"] = Intensity_dict[eachGroup][eachGroup+"_Intensity"] * this_ratio
        #     if volcanoData.empty:
        #         volcanoData = Intensity_dict[eachGroup]
        #     else:
        #         print(volcanoData)
        #         volcanoData = pd.merge(volcanoData,Intensity_dict[eachGroup],how="inner")

        
        fig, box = self.plot_HYE_colored(
            Intensity_dict,
            plot_options=plot_options,
            saved_settings=saved_settings,
            username=username,
            data_program=data_object["run_metadata"]["Processing App"][0]
        )
        CSV_link = None
        SVG_link = None
        # if WRITE_OUTPUT:
        #     # create the file for donwnload
        #     img_dir = os.path.join(APPFOLDER, "images/")
        #     if not os.path.exists(img_dir):
        #         Path(img_dir).mkdir(parents=True)

        #     fig.write_image(os.path.join(
        #         img_dir, f"{username}_abundance_volcano_Plot.svg"), format = "svg", validate = False, engine = "kaleido")
        #     # create the download CSV and its link

        #     data_dir = os.path.join(APPFOLDER, "csv/")
        #     if not os.path.exists(data_dir):
        #         Path(data_dir).mkdir(parents=True)
        #     volcanoData.to_csv(os.path.join(
        #         data_dir, f"{username}_up_down_regulated_volcano.csv"),
        #         index=False)
        #     print("Downloading links...")
        #     CSV_link = f"/files/{url_base}/csv/" \
        #         f"{username}_up_down_regulated_volcano.csv"

        #     # download SVG link
        #     SVG_link = f"/files/{url_base}/images/" \
        #         f"{username}_abundance_volcano_Plot.svg"
        
        return fig, box, CSV_link, SVG_link


    def plot_HYE_colored(self,allData_dict,
                            plot_options=None, 
                            saved_settings=None,
                            username=None,
                            data_program=None):
        
        fig = px.scatter(
            width=plot_options["width"],
            height=plot_options["height"],)
        box = px.box(
            width=plot_options["width"],
            height=plot_options["height"],
        )
        
        i = 0
        for each_organism in saved_settings["Order@"+plot_options["Group_By_Color"]]:
            groups = saved_settings["Order@"+plot_options["Group_By_Y"]]
            top = groups[0]
            bottom = groups[1]
            label = top+"/"+bottom
        
            these_groups=[]
            for each_group in  saved_settings.keys():
                if ("Order@" not in each_group and
                    saved_settings[each_group][plot_options["Group_By_Color"]] == each_organism and 
                    (saved_settings[each_group][plot_options["Group_By_Y"]] == top or saved_settings[each_group][plot_options["Group_By_Y"]] == bottom)):
                    these_groups.append(each_group)
                    if saved_settings[each_group][plot_options["Group_By_Y"]]==top:
                        top_group = each_group
                        top_intensity = each_group+"_Intensity"
                    elif saved_settings[each_group][plot_options["Group_By_Y"]]==bottom:
                        bottom_group = each_group
                        bottom_intensity = each_group+"_Intensity"
            
            x_axis = bottom_intensity

            this_data = pd.DataFrame()
            for each_group in these_groups:
                if this_data.empty:
                    this_data = allData_dict[each_group]
                else:
                    this_data = pd.merge(this_data,allData_dict[each_group],how="inner")
            this_data["organism"] = each_organism
            this_data["data_program"] = data_program

            if this_data.shape[0] != 0:
                print(each_organism)
                print(np.nanmedian(this_data[top_intensity]-this_data[bottom_intensity]),)
                fig.add_scatter(x=this_data[x_axis],
                                y=this_data[top_intensity]-this_data[bottom_intensity],
                                text=this_data["Accession"],
                                mode="markers", marker=dict(
                                    color=plot_options["colors"][i]),
                                name = each_organism)
                fig.add_hline(y=np.log2(saved_settings[top_group][plot_options["Group_By_Amount"]]/
                                    saved_settings[bottom_group][plot_options["Group_By_Amount"]]))
                box.add_box(x=this_data["data_program"],
                            y=this_data[top_intensity]-this_data[bottom_intensity],
                            fillcolor=plot_options["colors"][i],
                            name = each_organism)
                box.add_hline(y=np.log2(saved_settings[top_group][plot_options["Group_By_Amount"]]/
                                    saved_settings[bottom_group][plot_options["Group_By_Amount"]]))
            i = i + 1
        
        fig.update_traces(
                mode="markers",
                hovertemplate="%{text}<br>x=: %{x}"
                " <br>y=: %{y}")
        if plot_options["title"] != "" or plot_options["title"] is not None:
            plot_title = plot_options["title"] + " " + label
        else:
            plot_title = None
        if not plot_options["xlimits"] or plot_options["xlimits"] == "[]" or \
                not isinstance(plot_options["xlimits"], list):
            xlimits = None
        else:
            xlimits = plot_options["xlimits"]

        if not plot_options["ylimits"] or plot_options["ylimits"] == "[]" or \
                not isinstance(plot_options["ylimits"], list):
            ylimits = None
        else:
            ylimits = plot_options["ylimits"]

        fig.update_layout(
            font=plot_options["font"],

            title=plot_title,
            xaxis=dict(title=dict(
                text=plot_options["X Title"]), range=xlimits),
            yaxis=dict(title=dict(
                text=plot_options["Y Title"]), range=ylimits),
            plot_bgcolor='white',
            paper_bgcolor='white',

        )
        box.update_layout(
            font=plot_options["font"],

            title=plot_title,
            xaxis=dict(title=dict(
                text=plot_options["X Title"]), range=xlimits),
            yaxis=dict(title=dict(
                text=plot_options["Y Title"]), range=ylimits),
            plot_bgcolor='white',
            paper_bgcolor='white',

        )

        
        return fig, box
    # ###Ranked Abundance Plot####
    def Rank_Abundance_Plots(self, data_object,  plot_options, saved_settings, username=None):
        """_Prepare data for creating intensity volcano plots (two groups)_
        """
        group_names = []

        # Ranked Proteins
        ref_name = plot_options["reference_group"]

        run_ids = saved_settings[ref_name]["records"]
        if self.processor.data_type == "TMT":
            channel_ids = data_object["run_metadata"].loc[data_object["run_metadata"]["Run Identifier"].isin(run_ids)]["Channel Identifier"].to_list()
            ref_data = self.processor.filter_by_channel_id(data_object,channel_ids)["protein_abundance"] 
        elif self.processor.data_type == "LF":

            ref_data = self.processor.filter_by_id(data_object,list(run_ids))["protein_abundance"]
        ref_data["average_intensity"] = ref_data.mean(axis=1)
        ref_data = ref_data.sort_values(by="average_intensity",ascending = False).drop_duplicates(subset='Accession', keep='first').dropna().reset_index()
        ranked_proteins = pd.DataFrame({"Accession": ref_data["Accession"],
                                        "Rank": ref_data.index})

        # no compare groups is provided, compare first two
        for each_key in plot_options["compare groups"]:
            if each_key in plot_options["compare groups"] and saved_settings["Order@Conditions"]:
                group_names.append(each_key)

        # import the data
        group_dict = {}

        # filter runs into different groups
        i = 0
        runname_list = []  # contain list of run names list for each groups
        for eachGroup in group_names:
            run_ids = saved_settings[eachGroup]["records"]
            runname_sublist = data_object["run_metadata"].loc[data_object["run_metadata"]["Run Identifier"].isin(run_ids)]["Run Names"].to_list()
            if self.processor.data_type == "TMT":
                channel_ids = data_object["run_metadata"].loc[data_object["run_metadata"]["Run Identifier"].isin(run_ids)]["Channel Identifier"].to_list()
                group_dict[eachGroup] =  self.processor.filter_by_channel_id(
                    data_object,
                    channel_ids)  
            elif self.processor.data_type == "LF":

                group_dict[eachGroup] =  self.processor.filter_by_id(
                    data_object,
                    run_ids) # prevent the list from being changed
            runname_list.append(runname_sublist)
            i += 1
        # create a dictionary to store the intensity data
        Intensity_dict = {}

        for eachGroup in group_names:
            current_condition_data = self.processor.filter_by_missing_values(
                group_dict[eachGroup])
            current_prot_data = current_condition_data["protein_abundance"]
            current_prot_data["Average Reporter Ion Intensity"] = current_prot_data.mean(axis=1)
            current_prot_data = pd.merge(current_prot_data,ranked_proteins,how="inner")
            current_prot_data = current_prot_data.sort_values("Rank")
            Intensity_dict[eachGroup] = current_prot_data
            
            
        
        """ group1Data
                group1_Intensity  group1_stdev  group1_num   Accession
        0        2.824766e+05  1.708060e+05          15  A0A0B4J2D5
        1        2.650998e+06  6.259645e+05          15      A2RUR9
        2        1.973150e+05  5.645698e+04          15      A8MTJ3
        3        2.524020e+05  1.355699e+05          15      A8MWD9
        """


        fig = self.plot_ranked_abundance(Intensity_dict, plot_options, username)

        if self.write_output:
            # create the file for donwnload
            img_dir = os.path.join(self.app_folder, "images/")
            if not os.path.exists(img_dir):
                Path(img_dir).mkdir(parents=True)

            fig.write_image(os.path.join(
                img_dir, f"{username}_ranked_abundance_Plot.svg"), format = "svg", validate = False, engine = "kaleido")
            # create the download CSV and its link

            data_dir = os.path.join(self.app_folder, "csv/")
            if not os.path.exists(data_dir):
                Path(data_dir).mkdir(parents=True)
            for eachGroup in group_names:
                Intensity_dict[eachGroup].to_csv(os.path.join(data_dir, f"{username}_{eachGroup}_ranked_abundance.csv"),index=False)
            print("Downloading links...")
            CSV_link = f"/files/{self.url_base}/csv/" \
                f"{username}_up_down_regulated_volcano.csv"

            # download SVG link
            SVG_link = f"/files/{self.url_base}/images/" \
                f"{username}_abundance_volcano_Plot.svg"
            
            return fig, CSV_link, SVG_link


    def plot_ranked_abundance(self, allData,
                            plot_options=None,
                            username=None,):
        
        fig = px.scatter(
            width=plot_options["width"],
            height=plot_options["height"],)
        i = 0
        for eachCondition in allData:
            if allData[eachCondition].shape[0] != 0:
                fig.add_scatter(name=eachCondition,
                                x=allData[eachCondition]["Rank"],
                                y=allData[eachCondition]["Average Reporter Ion Intensity"],
                                text=allData[eachCondition]["Accession"],
                                mode="markers", 
                                marker_color = plot_options["color"][i])
            i = i + 1
        
        if plot_options["title"] != "" or plot_options["title"] is not None:
            plot_title = plot_options["title"]
        else:
            plot_title = None
        if not plot_options["xlimits"] or plot_options["xlimits"] == "[]" or \
                not isinstance(plot_options["xlimits"], list):
            xlimits = None
        else:
            xlimits = plot_options["xlimits"]

        if not plot_options["ylimits"] or plot_options["ylimits"] == "[]" or \
                not isinstance(plot_options["ylimits"], list):
            ylimits = None
        else:
            ylimits = plot_options["ylimits"]

        fig.update_layout(
            font=plot_options["font"],

            title=plot_title,
            xaxis=dict(title=dict(
                text=plot_options["X Title"]), range=xlimits),
            yaxis=dict(title=dict(
                text=plot_options["Y Title"]), range=ylimits),
            plot_bgcolor='white',
            paper_bgcolor='white',

        )

        
        return fig
    

    # ###Volcano plots####
    def volcano_plots(self,data_object,  plot_options, saved_settings, username=None, missing_values_max=33):
        """_Prepare data for creating intensity volcano plots (two groups)_
        """
        group_names = []

        # no compare groups is provided, compare first two
        for each_key in plot_options["compare groups"]:
            if each_key in plot_options["compare groups"] and saved_settings["Order@Conditions"]:
                group_names.append(each_key)

        # import the data
        group_dict = {}

        # filter runs into different groups
        i = 0
        runname_list = []  # contain list of run names list for each groups
        for eachGroup in group_names:
            run_ids = saved_settings[eachGroup]["records"]
            runname_sublist = data_object["run_metadata"].loc[data_object["run_metadata"]["Run Identifier"].isin(run_ids)]["Run Names"].to_list()
            if self.processor.data_type == "TMT":
                channel_ids = data_object["run_metadata"].loc[data_object["run_metadata"]["Run Identifier"].isin(run_ids)]["Channel Identifier"].to_list()
                group_dict[eachGroup] =  self.processor.filter_by_channel_id(
                    data_object,
                    channel_ids)  
            elif self.processor.data_type == "LF":

                group_dict[eachGroup] =  self.processor.filter_by_id(
                    data_object,
                    run_ids) # prevent the list from being changed
            runname_list.append(runname_sublist)
            i += 1
        # create a dictionary to store the intensity data
        Intensity_dict = {}

        for eachGroup in group_names:
            current_condition_data = self.processor.filter_by_missing_values(
                group_dict[eachGroup], missing_value_thresh=missing_values_max)
            Intensity_dict[eachGroup] =  current_condition_data["protein_abundance"]
            
        
        group1 = group_names[0]
        group2 = group_names[1]
        # calculate mean, standard deviation, and the number of non-null
        # elements for each row/protein
        group1Data = (Intensity_dict[group1]
                    .assign(**{group1+'_Intensity': Intensity_dict[group1].drop(
            columns=['Accession']).mean(axis=1),
            "group1_stdev":Intensity_dict[group1].drop(
                        columns=['Accession']).std(axis=1),
            "group1_num":Intensity_dict[group1].drop(
                        columns=['Accession']).shape[1] - Intensity_dict[
                        group1].isna().sum(axis=1)})
                    .loc[:, [group1+'_Intensity',
                            'group1_stdev',
                            'group1_num',
                            'Accession']])
        """ group1Data
                group1_Intensity  group1_stdev  group1_num   Accession
        0        2.824766e+05  1.708060e+05          15  A0A0B4J2D5
        1        2.650998e+06  6.259645e+05          15      A2RUR9
        2        1.973150e+05  5.645698e+04          15      A8MTJ3
        3        2.524020e+05  1.355699e+05          15      A8MWD9
        """
        group2Data = (Intensity_dict[group2]
                    .assign(**{group2+'_Intensity': Intensity_dict[group2].drop(
            columns=['Accession']).mean(axis=1),
            "group2_stdev":Intensity_dict[group2].drop(
                        columns=['Accession']).std(axis=1),
            "group2_num":Intensity_dict[group2].drop(
                        columns=['Accession']).shape[1] - Intensity_dict[
                        group2].isna().sum(axis=1)})
                    .loc[:, [group2+'_Intensity',
                            'group2_stdev',
                            'group2_num',
                            'Accession']])
        # find common proteins
        commonProts = group1Data[group1Data['Accession'].isin(group2Data['Accession'])].drop_duplicates().loc[:,"Accession"]
        print(commonProts.shape)
        print(group1Data.shape)
        print(group2Data.shape)


        # only leave common proteins
        group2Data = group2Data[group2Data['Accession'].isin(commonProts)].drop_duplicates()
        group1Data = group1Data[group1Data['Accession'].isin(commonProts)].drop_duplicates()
        print(group1Data.shape)
        print(group2Data.shape)

        group2Median = group2Data[group2+'_Intensity'].median(
            )
        group1Median = group1Data[group1+'_Intensity'].median(
            )
        #this is the median for each protein across both groups, you could do a ratio, but this puts it in terms of
        allmedian = pd.DataFrame({"col2":group2Data[group2+'_Intensity'],"col1":group1Data[group1+'_Intensity']}).median(axis=1,numeric_only=True)
        
        if (Intensity_dict[group1].shape[1] > 3 and
            Intensity_dict[group2].shape[1] > 3 and
                group2 != group1):
            # calculate the ratio between two group median,
            # will be used to normalize them
            ratio1 = allmedian / group1Median
            ratio2 = allmedian / group2Median

            # merge these two set of data together, adjust groups with ratio to 
            # median of all. Calculate pOriginal, p, significant
            # pOriginal is a numpy array or list of p-values
            # method is the method to be used for adjusting the p-values
            print(group2Data.columns)
            print(group1Data.columns)
            print(commonProts)
            volcanoData = (group2Data
                        .merge(group1Data, on='Accession', how='inner'))

            volcanoData = (volcanoData
                        .assign(**{group1+'_Intensity':lambda x: volcanoData[
                            group1+'_Intensity'] * ratio1}))
            volcanoData = (volcanoData
                        .assign(**{group2+'_Intensity':lambda x: volcanoData[
                            group2+'_Intensity'] * ratio2}))

            volcanoData = (volcanoData
                        .assign(
                            pOriginal=self.processor.t_test_from_summary_stats(
                                m1=volcanoData[group2+'_Intensity'],
                                m2=volcanoData[group1+'_Intensity'],
                                s1=volcanoData['group2_stdev'],
                                s2=volcanoData['group1_stdev'],
                                n1=volcanoData['group2_num'],
                                n2=volcanoData['group1_num'])))
            # filter out rows in volcanoData that have pOriginal == nan
            # if pOriginal is nan, then the p value will be nan
            volcanoData = volcanoData[volcanoData['pOriginal'].notna()]
            volcanoData = (volcanoData
                        .assign(benjamini=multipletests(volcanoData[
                            "pOriginal"], method='fdr_bh')[1]))

            volcanoData = (volcanoData
                        .assign(significant=volcanoData['benjamini'] < 0.01))

            # add upRegulated, downRegulated, and notRegulated columns
            volcanoData = volcanoData.assign(upRegulated=lambda x: (
                np.log10(volcanoData[group2+'_Intensity']/volcanoData[
                    group1+'_Intensity']) >np.log10(1.25)) & (volcanoData['significant']))

            volcanoData = volcanoData.assign(downRegulated=lambda x: (
                np.log10(volcanoData[group2+'_Intensity']/volcanoData[
                    group1+'_Intensity']) < -np.log10(1.25)) & (volcanoData['significant']))
            volcanoData = volcanoData.assign(notRegulated=lambda x: (abs(
            np.log10(volcanoData[group2+'_Intensity']/volcanoData[
                    group1+'_Intensity'])) <=np.log10(1.25)) & (~volcanoData['significant']))
            print(volcanoData)
            fig = self.plot_volcano_colored(
                volcanoData,
                label=f"({group2}/{group1})",
                saved_settings=saved_settings,
                plot_options=plot_options,
                username=username,
            )
            CSV_link = None
            SVG_link = None
            if self.write_output:
                # create the file for donwnload
                img_dir = os.path.join(self.app_folder, "images/")
                if not os.path.exists(img_dir):
                    Path(img_dir).mkdir(parents=True)

                fig.write_image(os.path.join(
                    img_dir, f"{username}_abundance_volcano_Plot.svg"), format = "svg", validate = False, engine = "kaleido")
                # create the download CSV and its link

                data_dir = os.path.join(self.app_folder, "csv/")
                if not os.path.exists(data_dir):
                    Path(data_dir).mkdir(parents=True)
                volcanoData.to_csv(os.path.join(
                    data_dir, f"{username}_up_down_regulated_volcano.csv"),
                    index=False)
                print("Downloading links...")
                CSV_link = f"/files/{self.url_base}/csv/" \
                    f"{username}_up_down_regulated_volcano.csv"

                # download SVG link
                SVG_link = f"/files/{self.url_base}/images/" \
                    f"{username}_abundance_volcano_Plot.svg"
            
            return fig, CSV_link, SVG_link


    def plot_volcano_colored(self, allData,
                            label, saved_settings,
                            plot_options=None,
                            username=None,):
        group_names = []

        # no compare groups is provided, compare first two
        for each_key in plot_options["compare groups"]:
            if each_key in plot_options["compare groups"] and saved_settings["Order@Conditions"]:
                group_names.append(each_key)
        
        total_labels = []
        left = group_names[0]+'_Intensity'
        right = group_names[1]+'_Intensity'
        downData = allData[allData['downRegulated']
                        == True]
        upData = allData[allData['upRegulated'] == True]

        fig = px.scatter(
            width=plot_options["width"],
            height=plot_options["height"],)
        if allData.shape[0] != 0:
            fig.add_scatter(x=np.log2(allData[right]/allData[left]),
                            y=-np.log10(allData["benjamini"]),
                            text=allData["Accession"],
                            mode="markers", marker=dict(
                                color=plot_options["all color"],size=20))
        if downData.shape[0] != 0:
            fig.add_scatter(x=np.log2(downData[right]/downData[left]),
                            y=-np.log10(downData["benjamini"]),
                            text=downData["Accession"],
                            mode="markers",
                            marker=dict(color=plot_options["down color"],size=20))
        if upData.shape[0] != 0:
            fig.add_scatter(x=np.log2(upData[right]/upData[left]),
                            y=-np.log10(upData["benjamini"]),
                            text=upData["Accession"],
                            mode="markers",
                            marker=dict(color=plot_options["up color"],size=20))
            fig.update_traces(
                mode="markers",
                hovertemplate="%{text}<br>x=: %{x}"
                " <br>y=: %{y}")
        fig.add_hline(y=2)
        fig.add_vline(x=-np.log2(10.0**-np.log10(1.25)))
        fig.add_vline(x=np.log2(10.0**-np.log10(1.25)))
        if plot_options["title"] != "" or plot_options["title"] is not None:
            plot_title = plot_options["title"] + " " + label
        else:
            plot_title = None
        if not plot_options["xlimits"] or plot_options["xlimits"] == "[]" or \
                not isinstance(plot_options["xlimits"], list):
            xlimits = None
        else:
            xlimits = plot_options["xlimits"]

        if not plot_options["ylimits"] or plot_options["ylimits"] == "[]" or \
                not isinstance(plot_options["ylimits"], list):
            ylimits = None
        else:
            ylimits = plot_options["ylimits"]

        fig.update_layout(
            font=plot_options["font"],

            showlegend=False,
            title=plot_title,
            xaxis=dict(title=dict(
                text=plot_options["X Title"]), range=xlimits),
            yaxis=dict(title=dict(
                text=plot_options["Y Title"]), range=ylimits),
            annotations=total_labels,
            plot_bgcolor='white',
            paper_bgcolor='white',

        )

        
        return fig


    # ###PCA plots####
    def PCA_plots(self, data_object, plot_options, saved_settings,username=None):
        """_Prepare data for creating intensity PCA plots (two groups)_
        """
        group_names = []

        # no compare groups is provided, compare first two
        for each_key in plot_options["compare groups"]:
            if each_key in plot_options["compare groups"] and saved_settings["Order@Conditions"]:
                group_names.append(each_key)
        # import the data
        group_dict = {}

        # filter runs into different groups
        i = 0
        runname_list = []  # this will contain list of run names list for each groups
        #print(saved_settings)
        for eachGroup in group_names:
            run_ids = saved_settings[eachGroup]["records"]
            runname_sublist = data_object["run_metadata"].loc[data_object["run_metadata"]["Run Identifier"].isin(run_ids)]["Run Names"].to_list()
            
            if self.processor.data_type == "TMT":
                channel_ids = data_object["run_metadata"].loc[data_object["run_metadata"]["Run Identifier"].isin(run_ids)]["Channel Identifier"].to_list()
                group_dict[eachGroup] =  self.processor.filter_by_channel_id(
                    data_object,
                    channel_ids)  
            elif self.processor.data_type == "LF":

                group_dict[eachGroup] =  self.processor.filter_by_id(
                    data_object,
                    run_ids)  # prevent the list from being changed
            runname_list.append(runname_sublist)
            i += 1       
            #print((runname_sublist))

        all_runs =[item for sublist in runname_list for item in sublist]

        # combined the data after filtering
        #  missing values and log2 transformation/normalization
        combined_infodata = pd.DataFrame() # store run names and group names
        combined_pcaData = pd.DataFrame() # store normalized data and protein names
        for eachGroup in group_names:

            current_condition_data = self.processor.filter_by_missing_values(
                group_dict[eachGroup])
            normalized_data = self.processor.NormalizeToMedian(
                current_condition_data["protein_abundance"],apply_log2=True)
            if self.data_type == "TMT":
                toFileDict = dict(zip(data_object["run_metadata"]["Channel Identifier"],data_object["run_metadata"]["Run Names"]))
            elif self.data_type == "LF":
                toFileDict = dict(zip(data_object["run_metadata"]["Run Identifier"],data_object["run_metadata"]["Run Names"]))
            toFileDict = self.processor.generate_column_to_name_mapping(normalized_data.columns, toFileDict)
            normalized_data.rename(columns = toFileDict,inplace=True)

            combined_infodata= pd.concat([combined_infodata, pd.DataFrame({
                "Sample_Groups": normalized_data
                .drop(
                    "Accession", axis=1).rename(columns = toFileDict).columns,
                "Type": eachGroup})])
            
            '''for x in range(len(list(combined_infodata["Sample_Groups"]))):
                print(list(combined_infodata["Sample_Groups"])[x])
            '''
            if combined_pcaData.empty:
                combined_pcaData = normalized_data
                print("Empty")
            else:
                combined_pcaData = pd.merge(combined_pcaData.drop_duplicates(), normalized_data.drop_duplicates())

        #normalize the data
        # using ratio of current group median value divide by the all groups median 
        # to create a scaling factor magicNUm to scale the each group
        quant_names = group_names
        while "Annotated Sequence" in quant_names:
            quant_names.remove("Annotated Sequence")
        while "Accession" in quant_names:
            quant_names.remove("Accession")

        while "Annotated Sequence" in all_runs:
            all_runs.remove("Annotated Sequence")
        while "Accession" in all_runs:
            all_runs.remove("Accession")

        # print(combined_pcaData.columns)
        # print(all_runs)
        # print(runname_list)
        # print(combined_pcaData[runname_list[0]].columns)
        # print(combined_pcaData[all_runs].columns)
        for n in range(len(quant_names)-2): #-2 because not Accession and Annotated Sequence 
            if "Annotated Sequence" in runname_list[n]:
                runname_list[n].remove("Annotated Sequence")
            if "Accession" in runname_list[n]:
                runname_list[n].remove("Accession")
            
            magicNum =np.nanmedian(combined_pcaData[runname_list[n]].dropna(how='all').to_numpy()
                                   ) /np.nanmedian(combined_pcaData[all_runs].dropna(how='all').to_numpy()) 
            for col in combined_pcaData[runname_list[n]].columns:
                combined_pcaData[col] = combined_pcaData[col]/magicNum


        #performs k-Nearest Neighbors imputation to fill in any missing values
        combined_pcaData = self.processor.impute_knn(combined_pcaData)
        combined_infodata.reset_index(drop=True, inplace=True)
        


        # perform PCA transform
        combined_pcaData, exp_var_pca = self.processor.CalculatePCA(combined_pcaData,
                                                        combined_infodata)

        return self.plot_PCA_plotly(combined_pcaData,
                            exp_var_pca,
                            plot_options=plot_options,
                            username=username,
                            )


    def plot_PCA_plotly(self, pca_panda,
                        exp_var_pca,
                        plot_options=None,
                        username=None,):

        CSV_link = None
        SVG_link = None

        # Assuming pca_data is a pandas dataframe containing PCA results
        # and "Type" is a column in the dataframe indicating the type of sample
        if not plot_options["xlimits"] or plot_options["xlimits"] == "[]" or \
                not isinstance(plot_options["xlimits"], list):
            xlimits = None
        else:
            xlimits = plot_options["xlimits"]

        if not plot_options["ylimits"] or plot_options["ylimits"] == "[]" or \
                not isinstance(plot_options["ylimits"], list):
            ylimits = None
        else:
            ylimits = plot_options["ylimits"]

        fig = px.scatter(pca_panda,
                        x='PC1',
                        y='PC2',
                        color="Type",
                        text="Sample_Groups",
                        symbol="Type",
                        color_discrete_sequence=plot_options["color"],

                        symbol_sequence=plot_options["symbol"],
                        size_max=30,
                        labels={'PC1': f'PC1 ({round(exp_var_pca[0]*100,2)}%)',
                                'PC2': f'PC2 ({round(exp_var_pca[1]*100,2)}%)',
                                'Type': 'Sample Type'}, title='PCA Plot',
                        width=plot_options["width"],
                        height=plot_options["height"],)

        fig.update_traces(
            mode="markers",
            marker=dict(size=plot_options["marker_size"],),
            hovertemplate="%{text}<br>PC1: %{x} <br>PC2: %{y}")
        fig.update_layout(
            plot_bgcolor="rgba(255, 255, 255, 255)",
            paper_bgcolor="rgba(255, 255, 255, 255)",
            font=plot_options["font"],
            title=plot_options["title"],
            xaxis=dict(linecolor='black',
                    showticklabels=False, mirror=True, range=xlimits),
            yaxis=dict(linecolor='black',
                    showticklabels=False, mirror=True, range=ylimits),
        )
        if self.write_output:
            # create the file for donwnload
            img_dir = os.path.join(self.app_folder, "images/")
            if not os.path.exists(img_dir):
                Path(img_dir).mkdir(parents=True)

            fig.write_image(os.path.join(
                img_dir, f"{username}_PCA_Plot.svg"), format = "svg", validate = False, engine = "kaleido")
            # create the download CSV and its link
            data_dir = os.path.join(self.app_folder, "csv/")
            if not os.path.exists(data_dir):
                Path(data_dir).mkdir(parents=True)
            pca_panda.to_csv(os.path.join(
                data_dir, f"{username}_PCA.csv"), index=False)
            print("Downloading links...")
            CSV_link = f"/files/{self.url_base}/csv/" \
                f"{username}_PCA.csv"

            # download SVG link
            SVG_link = f"/files/{self.url_base}/images/" \
                f"{username}_PCA_Plot.svg"

        return fig, CSV_link, SVG_link

    def heatmap_plots(self, data_object, plot_options, saved_settings, username=None):
        group_names = []

        # no compare groups is provided, compare first two
        for each_key in plot_options["compare groups"]:
            if each_key in plot_options["compare groups"] and saved_settings["Order@Conditions"]:
                group_names.append(each_key)
        # import the data
        group_dict = {}

        # filter runs into different groups
        i = 0
        runname_list = []  # this will contain list of run names list for each groups
        #print(saved_settings)
        for eachGroup in group_names:
            run_ids = saved_settings[eachGroup]["records"]
            runname_sublist = data_object["run_metadata"].loc[data_object["run_metadata"]["Run Identifier"].isin(run_ids)]["Run Names"].to_list()
            
            if self.processor.data_type == "TMT":
                channel_ids = data_object["run_metadata"].loc[data_object["run_metadata"]["Run Identifier"].isin(run_ids)]["Channel Identifier"].to_list()
                group_dict[eachGroup] =  self.processor.filter_by_channel_id(
                    data_object,
                    channel_ids)  
            elif self.processor.data_type == "LF":

                group_dict[eachGroup] =  self.processor.filter_by_id(
                    data_object,
                    run_ids)  # prevent the list from being changed
            runname_list.append(runname_sublist)
            i += 1       
            #print((runname_sublist))

        all_runs =[item for sublist in runname_list for item in sublist]

        # combined the data after filtering
        #  missing values and log2 transformation/normalization
        combined_infodata = pd.DataFrame() # store run names and group names
        combined_heatmap_data = pd.DataFrame() # store normalized data and protein names
        
        for eachGroup in group_names:

            current_condition_data = self.processor.filter_by_missing_values(
                group_dict[eachGroup])
            # print(current_condition_data)
            normalized_data = self.processor.NormalizeToMedian(
                current_condition_data["protein_abundance"],apply_log2=False) #apply this later
            if self.data_type == "TMT":
                toFileDict = dict(zip(data_object["run_metadata"]["Channel Identifier"],
                                [eachGroup + "_#" + str(i) for i in range(len(data_object["run_metadata"]["Channel Identifier"]))]))
            elif self.data_type == "LF":
                toFileDict = dict(zip(data_object["run_metadata"]["Run Identifier"],
                                [eachGroup + "_#" + str(i) for i in range(len(data_object["run_metadata"]["Run Identifier"]))]))
            print(toFileDict)
            toFileDict = self.processor.generate_column_to_name_mapping(normalized_data.columns, toFileDict)
            normalized_data.rename(columns = toFileDict,inplace=True)
            # print(normalized_data)

            combined_infodata= pd.concat([combined_infodata, pd.DataFrame({
                "Sample_Groups": normalized_data
                .drop(
                    "Accession", axis=1).rename(columns = toFileDict).columns,
                "Type": eachGroup})])
            
            '''for x in range(len(list(combined_infodata["Sample_Groups"]))):
                print(list(combined_infodata["Sample_Groups"])[x])
            '''
            if combined_heatmap_data.empty:
                combined_heatmap_data = normalized_data
                print("Empty")
            else:
                combined_heatmap_data = pd.merge(combined_heatmap_data, normalized_data)

        #normalize the data
        # using ratio of current group median value divide by the all groups median 
        # to create a scaling factor magicNUm to scale the each group
        quant_names = group_names
        while "Annotated Sequence" in quant_names:
            quant_names.remove("Annotated Sequence")
        while "Accession" in quant_names:
            quant_names.remove("Accession")

        while "Annotated Sequence" in all_runs:
            all_runs.remove("Annotated Sequence")
        while "Accession" in all_runs:
            all_runs.remove("Accession")

        for n in range(len(quant_names)-2): #-2 because not Accession and Annotated Sequence 
            if "Annotated Sequence" in runname_list[n]:
                runname_list[n].remove("Annotated Sequence")
            if "Accession" in runname_list[n]:
                runname_list[n].remove("Accession")
            
            magicNum =np.nanmedian(combined_heatmap_data[runname_list[
                n]].dropna(how='all').to_numpy()) /\
                    np.nanmedian(combined_heatmap_data[
                all_runs].dropna(how='all').to_numpy()) 
            for col in combined_heatmap_data[runname_list[
                    n]].columns:
                combined_heatmap_data[col] = combined_heatmap_data[col]/magicNum

            # print(combined_heatmap_data)
        if self.write_output == True or plot_options["significant_only"] == True:
            print(normalized_data)
            print(combined_heatmap_data)
            # How to transpose :(
            # detect differentially expressed proteins
            # reformat dataframe
            new_names = []
            for x in combined_heatmap_data.columns:
                if x.split(sep="_#")[0] in group_names:
                    new_names.append("Intensity" + x)
                else:
                    new_names.append(x)
            renamed_data = combined_heatmap_data.copy()
            renamed_data.columns = new_names
            # display(combined_heatmap_data)
            print(renamed_data)
            long_data = pd.wide_to_long(renamed_data.reset_index(),
                                        stubnames="Intensity",i="Accession",j="Sample",suffix=".*").reset_index()
            print(long_data)
            
            # if plot_options["agg_groups"]:
            #     # long_data= long_data.groupby(["Accession", "Group"]).agg({"Intensity": "mean"}).reset_index()

            #     # display(long_data)
            #     index = "Group"
            #     pivoted_data =  long_data.pivot(index=index,columns="Accession", values = "Intensity").reset_index()
            # else:
            #     index = "Sample"
            pivoted_data =  long_data.pivot(index="Sample",columns="Accession", values = "Intensity").reset_index()
            print(pivoted_data)
            pivoted_data["Group"] = pivoted_data["Sample"].str.replace("_#.*","",regex=True)
            # new_cols = combined_heatmap_data.columns.drop("Accession")
            # display(pivoted_data)
            pvalues = []
            means = []
            for col in pivoted_data.columns.drop(["Group","Sample"]):
                if (pivoted_data[col].dropna().shape[0]>=2):
                    pvalues.append(anova(*[pivoted_data.loc[pivoted_data["Group"]==x,col].dropna() for x in group_names])[0])
                    means.append([np.nanmean(pivoted_data.loc[pivoted_data["Group"]==x,col].dropna()) for x in group_names])
                    # print()
                else:
                    pvalues.append(np.nan)
                    means.append([np.nan])

            fold_changes = []
            # print(means)
            for eachprotein in means:
                length = len(eachprotein)
                total = np.nansum(eachprotein)
                current = []
                if type(eachprotein) != list:
                    print(str(eachprotein) + "******************")
                    current = [np.nan for x in group_names]
                else:
                    for eachMean in eachprotein:
                        each_fold_change = eachMean/((total-eachMean)/(length-1)) #don't include in comparison
                        current.append(each_fold_change)
                fold_changes.append(current)

            significances = []
            alpha = plot_options["alpha"] * 100 # anova function gives p values in %
            log_min_FC = np.log2(plot_options["min_fold_change"])
            for i in range(len(fold_changes)):
                current = False
                current_FCs = fold_changes[i]
                current_p = pvalues[i]
                if type(current_FCs) == list:
                    for each_FC in current_FCs:
                        if current_p < alpha and abs(np.log2(each_FC)) > log_min_FC:
                            current = True
                significances.append(current)
            combined_heatmap_data["adjusted_p_values"] = multipletests(pvalues, method='fdr_bh')[1] #adjust 
            combined_heatmap_data["fold_changes"] = fold_changes
            combined_heatmap_data["significant"] = significances



        if self.write_output:
            i = 0
            log_min_FC = np.log2(plot_options["min_fold_change"]) # anova function gives p values in %
            for eachGroup in group_names:
                if eachGroup not in combined_heatmap_data.columns:
                    j = 0
                    current = []
                    for x in combined_heatmap_data["fold_changes"]:
                        current.append(combined_heatmap_data["significant"][j] and abs(np.log2(x[i]))>log_min_FC)
                        j = j + 1
                    combined_heatmap_data[eachGroup] = current
                else:
                    print("ERROR: don't use group names that are also Accession numbers")
                i = i + 1

            data_dir = os.path.join(self.app_folder, "csv/")
            if not os.path.exists(data_dir):
                Path(data_dir).mkdir(parents=True)
            combined_heatmap_data.to_csv(os.path.join(
                data_dir, f"{username}_Heatmap.tsv"),sep="\t", index=False)
            
        if plot_options["significant_only"]:
            print(combined_heatmap_data.shape)
            combined_heatmap_data = combined_heatmap_data[combined_heatmap_data["significant"]]
            print(combined_heatmap_data.shape)
            combined_heatmap_data = combined_heatmap_data.drop(columns=["significant","adjusted_p_values","fold_changes"]+group_names)
        elif self.write_output:
            combined_heatmap_data = combined_heatmap_data.drop(columns=["significant","adjusted_p_values","fold_changes"]+group_names)

        if plot_options["log2_transform"]:
            for eachCol in combined_heatmap_data.columns.drop("Accession"):
                combined_heatmap_data[eachCol] = np.log2(combined_heatmap_data[eachCol])

        figure = px.imshow(combined_heatmap_data.set_index("Accession").dropna(), aspect="auto")
        CSV_link = None
        SVG_link = None
        return figure, CSV_link, SVG_link 