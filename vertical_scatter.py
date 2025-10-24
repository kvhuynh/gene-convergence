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
            "convergence_lengths": [],
            "convergence_count": None,
            "average_overlap_length": 0,
            "overlap_lengths": [],
            "overlap_count": 0
           }
    total_convergence = 0;
    total_overlap = 0;
    convergence_count = 0;
    overlap_count = 0;
    for element in data:
        converging_distance = data[element]["converging_distance"];
        if converging_distance > 0:
            total_convergence += converging_distance;
            res["convergence_lengths"].append(converging_distance);
        elif converging_distance <= 0:
            overlap_count += 1;
            total_overlap += converging_distance * -1;
            res["overlap_lengths"].append(converging_distance * -1)
        
        # an overlap is counted as convergent
        convergence_count += 1;
    # res["average_convergence_length"] = total_convergence / convergence_count;
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
        # === Vertical scatter plot for averages ===
        if plot_type == "con_avg":
            values = [
                VACV_convergence_data["convergence_lengths"],
                ASFV_convergence_data["convergence_lengths"],
                IIV_3_convergence_data["convergence_lengths"],
                IIV_6_convergence_data["convergence_lengths"],
                APMV_convergence_data["convergence_lengths"],
                MsV_convergence_data["convergence_lengths"]
            ]
        elif plot_type == "ovl_avg":
            values = [
                VACV_convergence_data["overlap_lengths"],
                ASFV_convergence_data["overlap_lengths"],
                IIV_3_convergence_data["overlap_lengths"],
                IIV_6_convergence_data["overlap_lengths"],
                APMV_convergence_data["overlap_lengths"],
                MsV_convergence_data["overlap_lengths"]
            ]

        colors = ['#2196F3', '#808080', '#9C27B0', '#F54927', '#FF9800', '#4CAF50']

        # Create vertical scatter with slight horizontal jitter for visibility
        for i, (vals, color) in enumerate(zip(values, colors)):
            x_jitter = np.random.normal(i, 0.05, size=len(vals))  # small spread in x
            plt.scatter(x_jitter, vals, alpha=0.7, color=color, s=20, label=labels[i])

            # Plot mean line
            mean_val = np.mean(vals)
            plt.hlines(mean_val, i - 0.2, i + 0.2, colors='black', linestyles='dashed', linewidth=1.5)
              # Add mean value as text above the line
            plt.text(i, mean_val + (0.02 * max(max(v) for v in values)), 
                     f"{mean_val:.2f}", ha='center', va='bottom', fontsize=9, fontweight='bold')

        plt.xticks(range(len(labels)), labels)
        plt.ylabel("Length")
        plt.title(f"{plot_type} (Vertical Scatter)")
        plt.ylim(0, max(max(v) for v in values) * 1.1)
        plt.legend(loc='upper right', frameon=False, fontsize=8)
        plt.tight_layout()
        plt.savefig(f"{plot_type}_20251024.png", dpi=300, bbox_inches='tight')
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
    plt.savefig(f"{plot_type}_20251024.png", dpi=300, bbox_inches='tight')
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






