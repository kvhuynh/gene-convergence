import matplotlib.pyplot as plt;
import numpy as np;
import json;

def parse_json(path):
    file = open(path).read();

    data = json.loads(file);
    return data;

def get_convergent_data(data):
    res = {
            "average_convergence_length": 0,
            "convergence_count": None,
            "average_overlap_length": 0,
            "overlap_count": 0
           }
    total_convergence = 0;
    total_overlap = 0;
    convergence_count = 0;
    overlap_count = 0;
    for element in data:
        converging_distance = data[element]["converging_distance"];
        if converging_distance > 0:
            convergence_count += 1;
            total_convergence += converging_distance;
        elif converging_distance <= 0:
            overlap_count += 1;
            total_overlap += converging_distance * -1;
    res["average_convergence_length"] = total_convergence / convergence_count;
    res["convergence_count"] = convergence_count; 

    res["average_overlap_length"] = total_overlap / overlap_count;
    res["overlap_count"] = overlap_count;
    return res;

def plot(MsV_convergence_data, VACV_convergence_data, IIV_3_convergence_data, IIV_6_convergence_data, APMV_convergence_data, ASFV_convergence_data, plot_type,
         MsV_total_genes=459, VACV_total_genes=218, IIV_3_total_genes=126, IIV_6_total_genes=463, APMV_total_genes=1018, ASFV_total_genes=168):

    plt.figure()
    # Desired left-to-right order:
    labels = ['VACV', 'ASFV', 'IIV-3', 'IIV-6', 'APMV', 'MsV']
    plt.title(f'{plot_type}')

    # Determine which values to plot
    if plot_type == "con_count":
        MsV_value = MsV_convergence_data["convergence_count"]
        VACV_value = VACV_convergence_data["convergence_count"]
        IIV_3_value = IIV_3_convergence_data["convergence_count"]
        IIV_6_value = IIV_6_convergence_data["convergence_count"]
        APMV_value = APMV_convergence_data["convergence_count"]
        ASFV_value = ASFV_convergence_data["convergence_count"]
    elif plot_type == "ovl_count":
        MsV_value = MsV_convergence_data["overlap_count"]
        VACV_value = VACV_convergence_data["overlap_count"]
        IIV_3_value = IIV_3_convergence_data["overlap_count"]
        IIV_6_value = IIV_6_convergence_data["overlap_count"]
        APMV_value = APMV_convergence_data["overlap_count"]
        ASFV_value = ASFV_convergence_data["overlap_count"]
    else:
        # For averages, keep your normal bar chart
        if plot_type == "con_avg":
            values = [
                VACV_convergence_data["average_convergence_length"],
                ASFV_convergence_data["average_convergence_length"],
                IIV_3_convergence_data["average_convergence_length"],
                IIV_6_convergence_data["average_convergence_length"],
                APMV_convergence_data["average_convergence_length"],
                MsV_convergence_data["average_convergence_length"]
            ]
        elif plot_type == "ovl_avg":
            values = [
                VACV_convergence_data["average_overlap_length"],
                ASFV_convergence_data["average_overlap_length"],
                IIV_3_convergence_data["average_overlap_length"],
                IIV_6_convergence_data["average_overlap_length"],
                APMV_convergence_data["average_overlap_length"],
                MsV_convergence_data["average_overlap_length"]
            ]

        plt.bar(labels, values, color=['#2196F3', '#808080', '#9C27B0', '#F54927', '#FF9800', '#4CAF50'])
        plt.ylim(0, max(values) * 1.15)
        for i, v in enumerate(values):
            plt.text(i, v + max(values) * 0.02, str(round(v, 2)),
                     ha='center', fontweight='bold')
        plt.tight_layout()
        plt.savefig(f"{plot_type}_20251023.png")
        plt.close()
        return

    # === For count plots (stacked bars) ===
    total_values = [VACV_total_genes, ASFV_total_genes, IIV_3_total_genes, IIV_6_total_genes, APMV_total_genes, MsV_total_genes]
    shaded_values = [VACV_value, ASFV_value, IIV_3_value, IIV_6_value, APMV_value, MsV_value]
    remaining_values = [t - s for t, s in zip(total_values, shaded_values)]

    plt.bar(labels, shaded_values,
            color=['#2196F3', '#808080', '#9C27B0', '#F54927', '#FF9800', '#4CAF50'],
            label='Actual Count')

    plt.bar(labels, remaining_values,
            bottom=shaded_values,
            color='lightgray',
            label='Remaining')

    for i, (val, total) in enumerate(zip(shaded_values, total_values)):
        plt.text(i, val / 2, str(val),
                 ha='center', va='center', color='white', fontweight='bold', fontsize=10)
        diff = total - val
        if diff > 0:
            plt.text(i, val + diff / 2, f"{diff}",
                     ha='center', va='center', color='black', fontsize=9)
                # Label total count above the bar
        plt.text(i, total + max(total_values) * 0.02,  # small offset above bar
                 f"{total}", ha='center', va='bottom', fontweight='bold', fontsize=10)


    plt.ylim(0, max(total_values) * 1.15)
    plt.ylabel('Gene Count', fontsize=11)
    plt.tight_layout()
    plt.savefig(f"{plot_type}_20251023.png", dpi=300, bbox_inches='tight')
    plt.close()



    # plt.savefig(f"{plot_type}.png") 
VACV = parse_json("./output/VACV_WR_convergence.json");
MsV = parse_json("./output/MsVT19_convergence.json");
IIV_3 = parse_json("./output/IIV-3_convergence.json");
IIV_6 = parse_json("./output/IIV-6_convergence.json");
APMV = parse_json("./output/APMV_convergence.json");
ASFV = parse_json("./output/ASFV_convergence.json");

MsV_convergence_data = get_convergent_data(MsV);
VACV_convergence_data = get_convergent_data(VACV);
IIV_3_convergence_data = get_convergent_data(IIV_3);
IIV_6_convergence_data = get_convergent_data(IIV_6);
APMV_convergence_data = get_convergent_data(APMV);
ASFV_convergence_data = get_convergent_data(ASFV);

plot(MsV_convergence_data, VACV_convergence_data, IIV_3_convergence_data, IIV_6_convergence_data, APMV_convergence_data, ASFV_convergence_data, "con_count");
plot(MsV_convergence_data, VACV_convergence_data, IIV_3_convergence_data, IIV_6_convergence_data,APMV_convergence_data, ASFV_convergence_data, "con_avg");
plot(MsV_convergence_data, VACV_convergence_data, IIV_3_convergence_data, IIV_6_convergence_data,APMV_convergence_data, ASFV_convergence_data, "ovl_count");
plot(MsV_convergence_data, VACV_convergence_data, IIV_3_convergence_data, IIV_6_convergence_data,APMV_convergence_data, ASFV_convergence_data, "ovl_avg");






