# -*- coding: utf-8 -*-
import os

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import yaml
from plotly.subplots import make_subplots

import utils

# Some information regarding the sankeys
# https://stackoverflow.com/questions/69498269/hide-plotly-sankey-nodes-and-links-while-preserving-total-node-value
# https://stackoverflow.com/questions/69463804/plotly-sankey-graph-data-formatting/69464558#69464558

# Colors from the tab10 color pallette
colors = {
    "Terrestrial solar irradiation": "rgba(180, 12, 13, opacity)",
    "Photosynthesis": "rgba(219, 97, 0, opacity)",
    "Biogas plant": "rgba(17, 128, 17, opacity)",
    "Reforming and conditioning": "rgba(116, 73, 156, opacity)",
    "Thermochemical gasification": "rgba(17, 128, 17, opacity)",
    "Conditioning": "rgba(116, 73, 156, opacity)",
    "Methanol_synth_conv": "rgba(5, 89, 147, opacity)",
    "Methanol_synth_direct": "rgba(5, 89, 147, opacity)",
    "loss": "rgba(97, 97, 97, opacity)",
    "PV": "rgba(219, 97, 0, opacity)",
    "EL+DAC": "rgba(154, 156, 7, opacity)",
}

opacity_node = str(0.95)
opacity_link = str(0.45)
COLOR_NODES = {k: v.replace("opacity", opacity_node) for k, v in colors.items()}
COLOR_LINKS = {k: v.replace("opacity", opacity_link) for k, v in colors.items()}
# get rid of extra color for loss-node
COLOR_NODES["loss"] = COLOR_LINKS["loss"]
STOP_VAL = 0.9

try:
    with open("./_variables.yml", "r") as vars:
        vars_dict = yaml.safe_load(vars)
except:
    with open("../_variables.yml", "r") as vars:
        vars_dict = yaml.safe_load(vars)


def remove_last_br(s):
    return s[::-1].replace(">rb<", " ", 1)[::-1]


def calc_flows_from_df(df):
    flows = []
    remainder_last_step = 100
    for remainder in df["Overall_efficiency_in_%"][1:]:
        loss = remainder_last_step - remainder
        flows.extend([loss, remainder])
        remainder_last_step = remainder
    flows.append(flows[-1])
    return flows


def add_legend(fig, label, color, y=0.0):
    # Add Legend:
    # Define custom legend items

    # Add custom legend as scatter plot points with transparent background
    for i, label in enumerate(label):
        fig.add_trace(
            go.Scatter(
                x=[None],  # Dummy x-coordinate
                y=[None],  # Dummy y-coordinate
                mode="markers",
                marker=dict(size=500, color=color[i], opacity=1),
                showlegend=True,
                name=label,
                hoverinfo="none",  # Disable hover info
            )
        )

    # Update layout for better visualization
    fig.update_layout(
        # title_text="Sankey Diagram with Custom Legend",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=y,
            xanchor="center",
            x=0.5,
            bgcolor="rgba(0,0,0,0)",  # Transparent background for the legend
        ),
    )


def layout_update(fig, fig_width=1800, fig_height=1000):
    fig.update_layout(
        autosize=False,
        width=fig_width,
        height=fig_height,
        margin=dict(t=50, l=10, r=65, b=0),
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        # showlegend=True,
        plot_bgcolor="rgba(0,0,0,0)",  # Transparent plot background
        font=dict(size=25),  # Node label font size
        legend=dict(font=dict(size=30)),  # Legend font size
    )
    fig.update_traces(textfont_color="black")


def layout_update_subplots(fig):
    fig.update_layout(
        annotations=[
            dict(
                font=dict(size=30),
                x=0,  # Position the text on the left
                xanchor="left",  # Aligns the text on the left
                y=1.012,  # Slightly above the plot area
                #     yanchor='bottom',     # Aligns the text at the bottom
                #     showarrow=False,      # No arrow for the annotation
                #     text="First Plot",    # Title text for the first plot
                #     xref='paper',         # x-reference is in paper coordinates
                #     yref='paper'          # y-reference is in paper coordinat
            ),
            dict(
                font=dict(size=30),
                x=0,  # Position the text on the left
                xanchor="left",  # Aligns the text on the left
                y=0.462,  # Slightly above the plot area
            ),
        ],  # Subplot title size
    )


def add_treemap(fig, domain=dict(x=[0.0, 0.12], y=[0.115, 1])):
    treemap_trace = go.Treemap(
        labels=["Solar energy", "Photosynthesis losses", " "],
        values=[100, 94.8, 5.2],
        parents=[
            "",
            "Solar energy",
            "Solar energy",
        ],
        marker_colors=[
            COLOR_NODES["Terrestrial solar irradiation"],
            COLOR_NODES["loss"],
            COLOR_NODES["Photosynthesis"],
        ],
        branchvalues="total",
        # textinfo="none",  # "percent root" # Hide the labels
        hoverinfo=None,
        textfont=dict(size=20),  # Increase the font size
        # Positioning domain of Treemap
        domain=domain,
    )

    # Add treemap trace to the figure
    fig.add_trace(treemap_trace)


def add_label_percentages(labels, flows, br=True):
    labs = []
    if br:
        br = "<br>"
    else:
        br = " "
    for idx, x in enumerate(labels):
        if x not in ["", "%"]:
            if idx == 0:
                labs.append(x + f"{br}100 %")
            elif idx % 2 == 1:
                labs.append(x + f"{br}{round(flows[idx+1], 1)} %")
            else:
                labs.append(x + f"{br}{round(flows[idx-1], 1)} %")
        elif x == "%":
            labs.append(f"{round(flows[idx+1], 1)} %")
        else:
            labs.append("")
    return labs


def create_x_vals(flows):
    x_vals = np.concatenate(
        [
            np.array([0.01]),
            np.repeat(
                np.linspace(0.01, STOP_VAL, num=(len(flows) // 2) + 1, endpoint=True),
                2,
            )[3:],
            np.array(
                [
                    0.99,
                ]
            ),
        ]
    )

    return x_vals


def create_y_vals(flows, start, stop):
    y_vals = [0.01] + [
        item
        for y_pos in np.linspace(start, stop, len(flows) // 2)
        for item in (0.01, y_pos)
    ]
    return y_vals


def plot_sankey_power_pathway(input_path, without_text=False):
    ## Sankey Power pathway ideal

    df_paths = [
        "PV_eff",
        "Gas_Provision_efficiencies_el_eff",
        "Methanol_synthesis_efficiencies_Electricity",
    ]

    dfs = []
    for df_path in df_paths:
        dfs.append(
            pd.read_csv(
                os.path.join(input_path, f"{df_path}.csv"), sep=";", index_col=0
            )
        )

    step_efficiencies = dfs[0]["Partial_efficiency_in_%"][1:].to_dict()
    step_efficiencies["Gas production (EL + DAC)"] = dfs[1].iloc[-1, -1]
    step_efficiencies["Methanol synthesis (direct)"] = dfs[2].iloc[-1, -1]
    step_efficiencies = {k: v / 100 for k, v in step_efficiencies.items()}

    df_electro = utils.table_creation("Electro", step_efficiencies)

    flows_electro = calc_flows_from_df(df_electro)
    flows_electro = [100] + flows_electro

    label_color_mapping = [
        ["Terrestrial solar<br>irradiation", "rgba(0,0,0,0)"],
        ["Absorbed photons", COLOR_NODES["Photosynthesis"]],
        ["Photons below band gap<br>not absorbed / transmitted.", COLOR_LINKS["loss"]],
        ["Energy", COLOR_NODES["Photosynthesis"]],
        ["Thermalization of photon<br>energy above band gap", COLOR_LINKS["loss"]],
        ["Electrical energy", COLOR_NODES["Photosynthesis"]],
        ["Extraction<br>losses", COLOR_LINKS["loss"]],
        ["Synthesis gas", COLOR_NODES["EL+DAC"]],
        [
            "Overpotentials in electrolysis<br> and energy for direct air capture",
            COLOR_LINKS["loss"],
        ],
        ["Methanol", COLOR_NODES["Methanol_synth_direct"]],
        ["Heat loss and<br>water separation", COLOR_LINKS["loss"]],
        ["", COLOR_LINKS["Methanol_synth_direct"]],
        ["", "rgba(0,0,0,0)"],
    ]

    if without_text:
        label_color_mapping = [["", x[1]] for x in label_color_mapping]

    labs = add_label_percentages([x[0] for x in label_color_mapping], flows_electro)
    x_vals = create_x_vals(flows_electro)
    y_vals = create_y_vals(flows_electro, 0.6, 0.15)
    y_vals[2] = 1.0
    y_vals[4] = 0.7
    y_vals[6] = 0.42
    y_vals[8] = 0.75
    y_vals[10] = 0.43

    df_nds = pd.DataFrame(
        {
            "labs": labs,
            "x_vals": x_vals,
            "y_vals": y_vals,
            "color": [x[1] for x in label_color_mapping],
        }
    )

    link_colors = [
        COLOR_LINKS["Terrestrial solar irradiation"],
        COLOR_LINKS["loss"],
        COLOR_LINKS["PV"],
        COLOR_LINKS["loss"],
        COLOR_LINKS["PV"],
        COLOR_LINKS["loss"],
        COLOR_LINKS["PV"],
        COLOR_LINKS["loss"],
        COLOR_LINKS["EL+DAC"],
        COLOR_LINKS["loss"],
        COLOR_LINKS["Methanol_synth_direct"],
        "rgba(0,0,0,0)",
    ]
    df_lnks = pd.DataFrame(
        {
            "source": np.concatenate(
                [
                    [0],
                    np.repeat(range(1, len(flows_electro) - 1, 2), 2),
                    [len(flows_electro) - 1],
                ]
            ),
            "target": list(range(1, len(flows_electro) + 1)),
            "value": flows_electro,
            "color": link_colors,
        }
    )

    sankey_1 = go.Sankey(
        arrangement="snap",
        # To get last label(s) on the right side:
        # Place on node right of the rest; have its line not show and everything else rgba(0,0,0,0). (Hover still exists however.)
        node={
            "label": df_nds["labs"],
            "x": df_nds["x_vals"],
            "y": df_nds["y_vals"],
            "color": df_nds["color"],
            "line": {"width": 0},
            "pad": 10,
        },  # 10 Pixels
        link={
            "source": df_lnks["source"],
            "target": df_lnks["target"],
            # Loss, remainder, ...
            # 100-partial eff, overall eff,
            "value": df_lnks["value"],
            "color": df_lnks["color"],
        },
    )

    ## Sankey Power pathway best_real

    step_efficiencies_electro_real = {
        "Power_extraction": vars_dict["eta_pv_real"] / 100,
        "Gas production (EL + DAC)": 0.56,
        "Methanol synthesis (direct)": vars_dict["eta_direct_meoh"] / 100,
    }

    df_electro_real = utils.table_creation("Electro", step_efficiencies_electro_real)

    flows_electro_real = calc_flows_from_df(df_electro_real)
    flows_electro_real = [100] + flows_electro_real

    label_color_mapping = [
        ["Terrestrial solar<br>irradiation", "rgba(0,0,0,0)"],
        ["Electrical energy", COLOR_NODES["Photosynthesis"]],
        [
            "Non absorption,<br>transmission,<br>thermalization,<br>and extraction<br>losses",
            COLOR_LINKS["loss"],
        ],
        ["Synthesis gas", COLOR_NODES["EL+DAC"]],
        [
            "Overpotentials in electrolysis<br> and energy for direct air capture",
            COLOR_LINKS["loss"],
        ],
        ["Methanol", COLOR_NODES["Methanol_synth_direct"]],
        ["Heat loss and<br>water separation", COLOR_LINKS["loss"]],
        ["", COLOR_LINKS["Methanol_synth_direct"]],
        ["", "rgba(0,0,0,0)"],
    ]

    if without_text:
        label_color_mapping = [["", x[1]] for x in label_color_mapping]

    labs = add_label_percentages(
        [x[0] for x in label_color_mapping], flows_electro_real
    )
    x_vals = list(df_nds["x_vals"])[:2] + list(df_nds["x_vals"])[6:]
    y_vals = create_y_vals(flows_electro_real, 0.6, 0.1)

    y_vals[2] = 0.85
    y_vals[4] = 0.63
    y_vals[6] = 0.33

    df_nds = pd.DataFrame(
        {
            "labs": labs,
            "x_vals": x_vals,
            "y_vals": y_vals,
            "color": [x[1] for x in label_color_mapping],
        }
    )

    link_colors = [
        COLOR_LINKS["Terrestrial solar irradiation"],
        COLOR_LINKS["loss"],
        COLOR_LINKS["PV"],
        COLOR_LINKS["loss"],
        COLOR_LINKS["EL+DAC"],
        COLOR_LINKS["loss"],
        COLOR_LINKS["Methanol_synth_direct"],
        "rgba(0,0,0,0)",
    ]
    df_lnks = pd.DataFrame(
        {
            "source": np.concatenate(
                [
                    [0],
                    np.repeat(range(1, len(flows_electro_real) - 1, 2), 2),
                    [len(flows_electro_real) - 1],
                ]
            ),
            "target": list(range(1, len(flows_electro_real) + 1)),
            "value": flows_electro_real,
            "color": link_colors,
        }
    )

    sankey_2 = go.Sankey(
        arrangement="snap",
        node={
            "label": df_nds["labs"],
            "x": df_nds["x_vals"],
            "y": df_nds["y_vals"],
            "color": df_nds["color"],
            "line": {"width": 0},
            "pad": 10,
        },  # 10 Pixels
        link={
            "source": df_lnks["source"],
            "target": df_lnks["target"],
            # Loss, remainder, ...
            # 100-partial eff, overall eff,
            "value": df_lnks["value"],
            "color": df_lnks["color"],
        },
    )

    # Create a subplot with 2 rows and 1 column of type 'domain'
    fig = make_subplots(
        rows=2,
        cols=1,
        specs=[[{"type": "domain"}], [{"type": "domain"}]],
        subplot_titles=("a) Theoretical efficiencies", "b) Present efficiencies"),
        vertical_spacing=0.16,
    )

    # Add the first Sankey diagram to the first subplot
    fig.add_trace(sankey_1, row=1, col=1)

    # Add the second Sankey diagram to the second subplot
    fig.add_trace(sankey_2, row=2, col=1)

    add_legend(
        fig,
        [
            "Solar cell",
            "Electrolysis & direct air capture",
            "Direct methanol synthesis",
            "Losses",
        ],
        [
            COLOR_NODES["PV"],
            COLOR_NODES["EL+DAC"],
            COLOR_NODES["Methanol_synth_direct"],
            COLOR_NODES["loss"],
        ],
        y=-0.16,
    )

    layout_update(fig)
    layout_update_subplots(fig)

    return fig


def plot_sankey_photosynthesis(input_path, without_text=False):
    # Sankey ideal photosynthesis

    # Read values flow values from csv
    df = pd.read_csv(
        os.path.join(input_path, "photosynthesis_eff.csv"), sep=";", index_col=0
    )

    flows_photosys = calc_flows_from_df(df)
    flows_photosys = [100] + flows_photosys

    label_color_mapping = [
        ["Terrestrial solar<br>irradiation", "rgba(0,0,0,0)"],
        ["%", COLOR_NODES["Photosynthesis"]],
        [
            "Energy outside the<br>photosynthetically<br>active spectral<br>range",
            COLOR_LINKS["loss"],
        ],
        ["%", COLOR_NODES["Photosynthesis"]],
        ["Reflection and<br>transmission", COLOR_LINKS["loss"]],
        ["%", COLOR_NODES["Photosynthesis"]],
        ["Photochemical<br>losses", COLOR_LINKS["loss"]],
        ["%", COLOR_NODES["Photosynthesis"]],
        ["Losses during<br>carbohydrate<br>synthesis", COLOR_LINKS["loss"]],
        ["%", COLOR_NODES["Photosynthesis"]],
        ["Photorespiratory<br>losses", COLOR_LINKS["loss"]],
        ["Energy stored<br>in biomass", COLOR_NODES["Photosynthesis"]],
        ["Respiration<br>losses", COLOR_LINKS["loss"]],
        ["", COLOR_LINKS["Photosynthesis"]],
        ["", "rgba(0,0,0,0)"],
    ]

    if without_text:
        label_color_mapping = [["", x[1]] for x in label_color_mapping]

    labs = add_label_percentages(
        [x[0] for x in label_color_mapping], flows_photosys, br=False
    )
    x_vals = create_x_vals(flows_photosys)
    y_vals = create_y_vals(flows_photosys, 0.6, 0.1)
    y_vals[2] = 0.83
    y_vals[4] = 0.65
    y_vals[6] = 0.55
    y_vals[8] = 0.39
    y_vals[10] = 0.28
    y_vals[12] = 0.20

    df_nds = pd.DataFrame(
        {
            "labs": labs,
            "x_vals": x_vals,
            "y_vals": y_vals,
            "color": [x[1] for x in label_color_mapping],
        }
    )

    link_colors = [
        COLOR_LINKS["Terrestrial solar irradiation"],
        COLOR_LINKS["loss"],
        COLOR_LINKS["PV"],
        COLOR_LINKS["loss"],
        COLOR_LINKS["PV"],
        COLOR_LINKS["loss"],
        COLOR_LINKS["PV"],
        COLOR_LINKS["loss"],
        COLOR_LINKS["PV"],
        COLOR_LINKS["loss"],
        COLOR_LINKS["PV"],
        COLOR_LINKS["loss"],
        COLOR_LINKS["PV"],
        "rgba(0,0,0,0)",
    ]
    df_lnks = pd.DataFrame(
        {
            "source": np.concatenate(
                [
                    [0],
                    np.repeat(range(1, len(flows_photosys) - 1, 2), 2),
                    [len(flows_photosys) - 1],
                ]
            ),
            "target": list(range(1, len(flows_photosys) + 1)),
            "value": flows_photosys,
            "color": link_colors,
        }
    )

    fig = go.Figure(
        go.Sankey(
            arrangement="snap",
            # To get last label(s) on the right side:
            # Place on node right of the rest; have its line not show and everything else rgba(0,0,0,0). (Hover still exists however.)
            node={
                "label": df_nds["labs"],
                "x": df_nds["x_vals"],
                "y": df_nds["y_vals"],
                "color": df_nds["color"],
                "line": {"width": 0},
                "pad": 10,
            },  # 10 Pixels
            link={
                "source": df_lnks["source"],
                "target": df_lnks["target"],
                # Loss, remainder, ...
                # 100-partial eff, overall eff,
                "value": df_lnks["value"],
                "color": df_lnks["color"],
            },
        )
    )

    layout_update(fig, fig_width=1800, fig_height=500)

    add_legend(
        fig,
        ["Photosynthesis", "Losses"],
        [
            COLOR_NODES["Photosynthesis"],
            COLOR_NODES["loss"],
        ],
        y=-0.1,
    )

    return fig


def add_sankey(
    fig,
    row,
    col,
    domain=dict(x=[0.0, 0.2], y=[0.115, 1]),
    x0=[0.115, 0.115],
    y0=[0.992, 0.96],
    x1=[0.145, 0.145],
    y1=[0.995, 0.08],
    remainder=5.18,
    line_width=4,
    without_text=False,
):
    labels = [
        "Terrestrial<br>solar<br>irradiation<br>(100)<br><br>",
        f"<br><br><br>Photo-<br>synthesis<br>loss<br>({100-remainder:.2f})",
        "",
        "",
    ]

    if without_text:
        labels = ["" for x in labels]

    extra_sankey_trace = go.Sankey(
        node={
            "label": labels,
            "x": [0, 0.1, 0.2, 0.2],
            "y": [0.01, 0.48, 0.03, 0.63],
            "color": [
                "rgba(0,0,0,0)",
                COLOR_NODES["Photosynthesis"],
                COLOR_NODES["Biogas plant"],
                COLOR_NODES["loss"],
            ],
            "line": dict(color="black", width=0),
            "pad": 15,
        },
        link={
            "source": [0, 1, 1],
            "target": [1, 2, 3],
            "value": [100, remainder, 100 - remainder],
            "color": [
                COLOR_LINKS["Terrestrial solar irradiation"],
                COLOR_LINKS["Photosynthesis"],
                COLOR_LINKS["loss"],
            ],
        },
        domain=domain,
    )

    # Add treemap trace to the figure
    fig.add_trace(extra_sankey_trace, row=row, col=col)


def plot_sankey_anaerobic_digestion(input_path, without_text=False):
    ## Sankey Biogas (anaerobic digestion) pathway ideal

    df_biogas = [
        "photosynthesis_eff",
        "AD_eff",
        "Biogas_Ref_efficiencies",
        "Methanol_synthesis_efficiencies_BIO",
    ]

    dfs_biogas = []
    for df_path in df_biogas:
        dfs_biogas.append(
            pd.read_csv(
                os.path.join(input_path, f"{df_path}.csv"), sep=";", index_col=0
            )
        )

    step_efficiencies = {}
    step_efficiencies["Photosynthesis"] = dfs_biogas[0].iloc[-1, -1]
    step_efficiencies["Biogas plant"] = dfs_biogas[1].iloc[-1, -1]
    step_efficiencies["Reforming and conditioning"] = dfs_biogas[2].iloc[-1, -1]
    step_efficiencies["Methanol synthesis (conventional)"] = dfs_biogas[3].iloc[-1, -1]

    step_efficiencies = {k: v / 100 for k, v in step_efficiencies.items()}

    df_bio_biogas = utils.table_creation("Biogas", step_efficiencies)

    flows_bio_biogas = calc_flows_from_df(df_bio_biogas)
    flows_bio_biogas = flows_bio_biogas[2:]

    labels = [
        "Energy stored<br>in biomass",
        "Heat loss and<br>microorganism reproduction",
        "",
        "CO<sub>2</sub> separation and O<sub>2</sub> provision",
        "",
        "Heat loss",
        "Energy stored<br>in methanol",
        "",
    ]

    labs = [
        (
            x
            + (
                f"<br>(5.18)"
                if idx == 0
                else f"<br>({round(flows_bio_biogas[idx-1], 2):.2f})"
            )
            if x
            else f"({round(flows_bio_biogas[idx-1], 2):.2f})"
        )
        for idx, x in enumerate(labels[:-1])
    ]
    # labs = ["" for x in labs]

    x_vals = np.append(
        np.repeat(
            np.linspace(
                0.25, STOP_VAL, num=(len(flows_bio_biogas) // 2) + 1, endpoint=True
            ),
            2,
        )[1:],
        0.99,
    )
    y_vals = [
        item
        for y_pos in np.linspace(0.6, 0.1, len(flows_bio_biogas) // 2 + 1)
        for item in (0.01, y_pos)
    ]
    # y_vals =[0.01, np.float64(0.6), 0.01, 0.9, 0.01, 0.32, 0.01, np.float64(0.1)]
    y_vals[1] = 0.99
    y_vals[3] = 0.9
    y_vals[5] = 0.73

    sankey_1 = go.Sankey(
        arrangement="snap",
        # To get last label(s) on the right side:
        # Place on node right of the rest; have its line not show and everything else white. (Hover still exists however.)
        node={
            "label": [""] * len(x_vals) if without_text else labs,
            "x": x_vals,
            "y": y_vals,
            "color": [
                COLOR_NODES["Biogas plant"],
                COLOR_NODES["loss"],
                COLOR_NODES["Reforming and conditioning"],
                COLOR_NODES["loss"],
                COLOR_NODES["Methanol_synth_conv"],
                COLOR_NODES["loss"],
                COLOR_LINKS["Methanol_synth_conv"],
                "rgba(0,0,0,0)",
            ],
            "line": dict(color="black", width=0),
            "pad": 15,
        },  # 10 Pixels
        link={
            "source": [
                i for x in range(0, len(flows_bio_biogas) - 1, 2) for i in (x, x)
            ]
            + [len(flows_bio_biogas) - 1],
            "target": [x for x in range(1, len(flows_bio_biogas) + 1)],
            # Loss, remainder, ...
            # 100-partial eff, overall eff,
            "value": flows_bio_biogas,
            "color": [
                COLOR_LINKS["loss"],
                COLOR_LINKS["Biogas plant"],
                COLOR_LINKS["loss"],
                COLOR_LINKS["Reforming and conditioning"],
                COLOR_LINKS["loss"],
                COLOR_LINKS["Methanol_synth_conv"],
                "rgba(0,0,0,0)",
            ],
        },
    )

    # ------------------------------------------------------

    ## Sankey Biogas (anaerobic digestion) pathway real best

    step_efficiencies = {
        "Photosynthesis": vars_dict["eta_photosynthesis_real"] / 100,
        "Biogas plant": vars_dict["eta_anaerobic_digestion_higher_value"] / 100,
        "Reforming and conditioning": vars_dict["eta_ref_conditioning_higher_value"]
        / 100,
        "Methanol synthesis (conventional)": vars_dict["eta_conv_meoh"] / 100,
    }

    df_bio_biogas = utils.table_creation("Biogas", step_efficiencies)

    flows_bio_biogas = calc_flows_from_df(df_bio_biogas)
    flows_bio_biogas = flows_bio_biogas[2:]

    labels = [
        "Energy stored<br>in biomass",
        "Heat loss and<br>microorganism reproduction",
        "",
        "CO<sub>2</sub> separation and O<sub>2</sub> provision",
        "",
        "Heat loss",
        "Energy stored<br>in methanol",
        "",
    ]

    labs = [
        (
            x
            + (
                f"<br>(1.00)"
                if idx == 0
                else f"<br>({round(flows_bio_biogas[idx-1], 2):.2f})"
            )
            if x
            else f"({round(flows_bio_biogas[idx-1], 2):.2f})"
        )
        for idx, x in enumerate(labels[:-1])
    ]
    # labs = ["" for x in labs]

    x_vals = np.append(
        np.repeat(
            np.linspace(
                0.25, STOP_VAL, num=(len(flows_bio_biogas) // 2) + 1, endpoint=True
            ),
            2,
        )[1:],
        0.99,
    )
    y_vals = [
        item
        for y_pos in np.linspace(0.6, 0.1, len(flows_bio_biogas) // 2 + 1)
        for item in (0.01, y_pos)
    ]
    # y_vals =[0.01, np.float64(0.6), 0.01, 0.9, 0.01, 0.32, 0.01, np.float64(0.1)]
    y_vals[1] = 0.8
    y_vals[3] = 0.65
    y_vals[5] = 0.40

    # Introduce not visible node for scaling effect.
    labs = ["", "", *labs]  #
    x_vals = [0, 0.15, *x_vals]  #
    y_vals = [0.01, 0.9, *y_vals]  #

    flows_bio_biogas = [4.2, 1.0, *flows_bio_biogas]
    source = [i for x in range(0, len(flows_bio_biogas) - 1, 2) for i in (x, x)] + [
        len(flows_bio_biogas) - 1
    ]
    target = [x for x in range(1, len(flows_bio_biogas) + 1)]

    sankey_2 = go.Sankey(
        arrangement="snap",
        # To get last label(s) on the right side:
        # Place on node right of the rest; have its line not show and everything else white. (Hover still exists however.)
        node={
            "label": [""] * len(x_vals) if without_text else labs,
            "x": x_vals,
            "y": y_vals,
            "color": [
                "rgba(0,0,0,0)",
                "rgba(0,0,0,0)",
                COLOR_NODES["Biogas plant"],
                COLOR_NODES["loss"],
                COLOR_NODES["Reforming and conditioning"],
                COLOR_NODES["loss"],
                COLOR_NODES["Methanol_synth_conv"],
                COLOR_NODES["loss"],
                COLOR_LINKS["Methanol_synth_conv"],
                "rgba(0,0,0,0)",
            ],
            "line": dict(color="black", width=0),
            "pad": 15,
        },  # 10 Pixels
        link={
            "source": source,
            "target": target,
            # Loss, remainder, ...
            # 100-partial eff, overall eff,
            "value": flows_bio_biogas,
            "color": [
                "rgba(0,0,0,0)",
                "rgba(0,0,0,0)",
                COLOR_LINKS["loss"],
                COLOR_LINKS["Biogas plant"],
                COLOR_LINKS["loss"],
                COLOR_LINKS["Reforming and conditioning"],
                COLOR_LINKS["loss"],
                COLOR_LINKS["Methanol_synth_conv"],
                "rgba(0,0,0,0)",
            ],
        },
    )

    # Create a subplot with 2 rows and 1 column of type 'domain'
    fig = make_subplots(
        rows=2,
        cols=1,
        specs=[[{"type": "domain"}], [{"type": "domain"}]],
        subplot_titles=("a) Theoretical efficiencies", "b) Present efficiencies"),
        vertical_spacing=0.13,
    )

    # Add the first Sankey diagram to the first subplot
    fig.add_trace(sankey_1, row=1, col=1)

    add_sankey(
        fig,
        row=1,
        col=1,
        domain=dict(x=[0.0, 0.11], y=[0.072, 1]),
        x0=[0.125, 0.125],
        y0=[0.995, 0.98],  # y-vals on left side
        y1=[0.998, 0.62],  # y vales on right side
        remainder=5.2,
        without_text=without_text,
    )

    # Add the second Sankey diagram to the second subplot
    fig.add_trace(sankey_2, row=2, col=1)

    add_sankey(
        fig,
        row=2,
        col=1,
        domain=dict(x=[0.0, 0.11], y=[0.072, 1]),
        x0=[0.125, 0.125],
        y0=[0.422, 0.421],  # y-vals on left side
        y1=[0.433, 0.368],  # y vales on right side
        remainder=1.0,
        without_text=without_text,
    )

    add_legend(
        fig,
        [
            "Photosynthesis",
            "Anaerobic digestion",
            "Reforming & conditioning",
            "Conventional methanol synthesis",
            "Losses",
        ],
        [
            COLOR_NODES["Photosynthesis"],
            COLOR_NODES["Biogas plant"],
            COLOR_NODES["Reforming and conditioning"],
            COLOR_NODES["Methanol_synth_conv"],
            COLOR_NODES["loss"],
        ],
        y=-0.11,
    )

    layout_update(fig)
    layout_update_subplots(fig)

    return fig


def plot_sankey_biomass_gasification(input_path, without_text):
    ## Sankey Biomass Reforming pathway ideal
    # -------------------------------------------------

    df_biogas = [
        "photosynthesis_eff",
        "AD_eff",
        "Biogas_Ref_efficiencies",
        "Methanol_synthesis_efficiencies_BIO",
    ]

    dfs_biogas = []
    for df_path in df_biogas:
        dfs_biogas.append(
            pd.read_csv(
                os.path.join(input_path, f"{df_path}.csv"), sep=";", index_col=0
            )
        )

    df_gasification = [
        "photosynthesis_eff",
        "Gasification_eff",
        "Gas_conditioning_efficiencies",
        "Methanol_synthesis_efficiencies_BIO",
    ]

    dfs_gasification = []
    for df_path in df_gasification:
        dfs_gasification.append(
            pd.read_csv(
                os.path.join(input_path, f"{df_path}.csv"), sep=";", index_col=0
            )
        )

    step_efficiencies = {}
    step_efficiencies["Photosynthesis"] = dfs_biogas[0].iloc[-1, -1]
    step_efficiencies["Thermochemical gasification"] = dfs_gasification[1].iloc[-1, -1]
    step_efficiencies["Conditioning"] = dfs_gasification[2].iloc[-1, -1]
    step_efficiencies["Methanol synthesis (conventional)"] = dfs_gasification[3].iloc[
        -1, -1
    ]
    step_efficiencies = {k: v / 100 for k, v in step_efficiencies.items()}

    df_bio_gasification = utils.table_creation("Gasification", step_efficiencies)

    flows_bio_gasification = calc_flows_from_df(df_bio_gasification)
    flows_bio_gasification = flows_bio_gasification[2:]

    labels = [
        "Energy stored<br>in biomass",
        "Thermochemical gasification",
        "",
        # "CO<sub>2</sub> separation and<br>O<sub>2</sub> provision",
        "CO<sub>2</sub> separation and O<sub>2</sub> provision",
        "",
        "Heat loss",
        "Energy stored<br>in methanol",
        "",
    ]
    labs = [
        (
            x
            + (
                f"<br>(5.18)"
                if idx == 0
                else f"<br>({round(flows_bio_gasification[idx-1], 2):.2f})"
            )
            if x
            else f"({round(flows_bio_gasification[idx-1], 2):.2f})"
        )
        for idx, x in enumerate(labels[:-1])
    ]
    labs[1] = remove_last_br(labs[1])
    labs[-4] = remove_last_br(labs[-4])

    if without_text:
        labs = ["" for x in labs]

    x_vals = np.append(
        np.repeat(
            np.linspace(
                0.25,
                STOP_VAL,
                num=(len(flows_bio_gasification) // 2) + 1,
                endpoint=True,
            ),
            2,
        )[1:],
        0.99,
    )
    y_vals = [
        item
        for y_pos in np.linspace(0.99, 0.1, len(flows_bio_gasification) // 2 + 1)
        for item in (0.01, y_pos)
    ]

    node_colors = [
        COLOR_NODES["Biogas plant"],
        COLOR_NODES["loss"],
        COLOR_NODES["Reforming and conditioning"],
        COLOR_NODES["loss"],
        COLOR_NODES["Methanol_synth_conv"],
        COLOR_NODES["loss"],
        COLOR_LINKS["Methanol_synth_conv"],
        "rgba(0,0,0,0)",
    ]

    source_nodes = [
        i for x in range(0, len(flows_bio_gasification) - 1, 2) for i in (x, x)
    ] + [len(flows_bio_gasification) - 1]
    target_nodes = [x for x in range(1, len(flows_bio_gasification) + 1)]

    flow_values = flows_bio_gasification

    link_colors = [
        COLOR_LINKS["loss"],
        COLOR_LINKS["Thermochemical gasification"],
        COLOR_LINKS["loss"],
        COLOR_LINKS["Conditioning"],
        COLOR_LINKS["loss"],
        COLOR_LINKS["Methanol_synth_conv"],
        "rgba(0,0,0,0)",
    ]

    # Getting rid of Thermochemical gasification loss:
    labs.pop(1)
    x_vals = np.delete(x_vals, 1)
    y_vals = np.delete(y_vals, 1)
    node_colors.pop(1)

    source_nodes = [0, 1, 1, 3, 3, 5]
    target_nodes.pop()
    flow_values.pop(0)
    link_colors.pop(0)

    y_vals[2] = 1

    sankey_1 = go.Sankey(
        arrangement="snap",
        # To get last label(s) on the right side:
        # Place on node right of the rest; have its line not show and everything else white. (Hover still exists however.)
        node={
            "label": labs,
            "x": x_vals,
            "y": y_vals,
            "color": node_colors,
            "line": dict(color="black", width=0),
            "pad": 15,
        },  # 10 Pixels
        link={
            "source": source_nodes,
            "target": target_nodes,
            # Loss, remainder, ...
            # 100-partial eff, overall eff,
            "value": flow_values,
            "color": link_colors,
        },
    )

    # ------------------------------------------------------

    ## Sankey Biomass Reforming pathway best real

    step_efficiencies = {
        "Photosynthesis": vars_dict["eta_photosynthesis_real"] / 100,
        "Thermochemical gasification": vars_dict["eta_gasification_higher_value"] / 100,
        "Conditioning": vars_dict["eta_gas_conditioning_higher_value"] / 100,
        "Methanol synthesis (conventional)": vars_dict["eta_conv_meoh"] / 100,
    }

    df_bio_gasification = utils.table_creation("Gasification", step_efficiencies)

    flows_bio_gasification = calc_flows_from_df(df_bio_gasification)
    flows_bio_gasification = flows_bio_gasification[2:]

    labels = [
        "Energy stored<br>in biomass",
        "Thermochemical gasification",
        "",
        # "CO<sub>2</sub> separation and<br>O<sub>2</sub> provision",
        "CO<sub>2</sub> separation and O<sub>2</sub> provision",
        "",
        "Heat loss",
        "Energy stored<br>in methanol",
        "",
    ]
    labs = [
        (
            x
            + (
                f"<br>(1.00)"
                if idx == 0
                else f"<br>({round(flows_bio_gasification[idx-1], 2):.2f})"
            )
            if x
            else f"({round(flows_bio_gasification[idx-1], 2):.2f})"
        )
        for idx, x in enumerate(labels[:-1])
    ]
    labs[1] = remove_last_br(labs[1])
    labs[-4] = remove_last_br(labs[-4])

    if without_text:
        labs = ["" for x in labs]

    x_vals = np.append(
        np.repeat(
            np.linspace(
                0.25,
                STOP_VAL,
                num=(len(flows_bio_gasification) // 2) + 1,
                endpoint=True,
            ),
            2,
        )[1:],
        0.99,
    )
    y_vals = [
        item
        for y_pos in np.linspace(0.65, 0.1, len(flows_bio_gasification) // 2 + 1)
        for item in (0.01, y_pos)
    ]

    node_colors = [
        "rgba(0,0,0,0)",
        "rgba(0,0,0,0)",
        COLOR_NODES["Biogas plant"],
        COLOR_NODES["loss"],
        COLOR_NODES["Reforming and conditioning"],
        COLOR_NODES["loss"],
        COLOR_NODES["Methanol_synth_conv"],
        COLOR_NODES["loss"],
        COLOR_LINKS["Methanol_synth_conv"],
        "rgba(0,0,0,0)",
    ]

    link_colors = [
        "rgba(0,0,0,0)",
        "rgba(0,0,0,0)",
        COLOR_LINKS["loss"],
        COLOR_LINKS["Thermochemical gasification"],
        COLOR_LINKS["loss"],
        COLOR_LINKS["Conditioning"],
        COLOR_LINKS["loss"],
        COLOR_LINKS["Methanol_synth_conv"],
        "rgba(0,0,0,0)",
    ]

    # Introduce not visible node for scaling effect.
    labs = ["", "", *labs]
    x_vals = [0, 0.15, *x_vals]
    y_vals = [0.01, 0.9, *y_vals]
    flows_bio_gasification = [4.2, 1.0, *flows_bio_gasification]
    source_nodes = [
        i for x in range(0, len(flows_bio_gasification) - 1, 2) for i in (x, x)
    ] + [len(flows_bio_gasification) - 1]
    target_nodes = [x for x in range(1, len(flows_bio_gasification) + 1)]

    sankey_2 = go.Sankey(
        arrangement="snap",
        # To get last label(s) on the right side:
        # Place on node right of the rest; have its line not show and everything else white. (Hover still exists however.)
        node={
            "label": labs,
            "x": x_vals,
            "y": y_vals,
            "color": node_colors,
            "line": dict(color="black", width=0),
            "pad": 15,
        },  # 10 Pixels
        link={
            "source": source_nodes,
            "target": target_nodes,
            # Loss, remainder, ...
            # 100-partial eff, overall eff,
            "value": flows_bio_gasification,
            "color": link_colors,
        },
    )

    # Create a subplot with 2 rows and 1 column of type 'domain'
    fig = make_subplots(
        rows=2,
        cols=1,
        specs=[[{"type": "domain"}], [{"type": "domain"}]],
        subplot_titles=("a) Theoretical efficiencies", "b) Present efficiencies"),
        vertical_spacing=0.13,
    )

    # Add the first Sankey diagram to the first subplot
    fig.add_trace(sankey_1, row=1, col=1)

    add_sankey(
        fig,
        row=1,
        col=1,
        domain=dict(x=[0.0, 0.11], y=[0.072, 1]),
        x0=[0.125, 0.125],
        y0=[0.995, 0.98],  # y-vals on left side
        y1=[0.998, 0.62],  # y vales on right side
        remainder=5.2,
        without_text=without_text,
    )

    # Add the second Sankey diagram to the second subplot
    fig.add_trace(sankey_2, row=2, col=1)

    add_sankey(
        fig,
        row=2,
        col=1,
        domain=dict(x=[0.0, 0.11], y=[0.072, 1]),
        x0=[0.125, 0.125],
        y0=[0.422, 0.421],  # y-vals on left side
        y1=[0.433, 0.368],  # y vales on right side
        remainder=1.0,
        without_text=without_text,
    )

    add_legend(
        fig,
        [
            "Photosynthesis",
            "Thermochemical gasification",
            "Conditioning",
            "Conventional methanol synthesis",
            "Losses",
        ],
        [
            COLOR_NODES["Photosynthesis"],
            COLOR_NODES["Biogas plant"],
            COLOR_NODES["Reforming and conditioning"],
            COLOR_NODES["Methanol_synth_conv"],
            COLOR_NODES["loss"],
        ],
        y=-0.14,
    )

    layout_update(fig)
    layout_update_subplots(fig)

    return fig


def main(input_path="./data/results/", output_path="./img/"):
    # Label positioning finetuned afterwards
    fig = plot_sankey_power_pathway(input_path=input_path, without_text=False)
    fig.write_image(os.path.join(output_path, "sankey_electro_w_text.png"))
    # fig.write_image(os.path.join(output_path,"sankey_electro_wo_text.png"))

    # Label positioning finetuned afterwards
    fig = plot_sankey_photosynthesis(input_path=input_path, without_text=False)
    fig.write_image(os.path.join(output_path, "sankey_photosynthesis_w_text.png"))
    # fig.write_image(os.path.join(output_path,"sankey_photosynthesis_wo_text.png"))

    # Label positioning finetuned afterwards
    fig = plot_sankey_anaerobic_digestion(input_path=input_path, without_text=False)
    fig.write_image(os.path.join(output_path, "sankey_digestion_w_text.png"))
    # fig.write_image(os.path.join(output_path,"sankey_digestion_wo_text.png"))

    # Label positioning finetuned afterwards
    fig = plot_sankey_biomass_gasification(input_path=input_path, without_text=False)
    fig.write_image(os.path.join(output_path, "sankey_gasification_w_text.png"))
    # fig.write_image(os.path.join(output_path,"sankey_gasification_wo_text.png"))


if __name__ == "__main__":
    main()
