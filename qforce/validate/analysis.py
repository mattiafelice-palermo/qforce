import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd


def plot_correlation(md_output, qm_output):
    common_configurations = list(set(qm_output.ids).intersection(md_output.ids))
    discarded_configurations = set(qm_output.ids) - set(md_output.ids)

    qm_energies = np.array([qm_output.spe[config - 1] for config in common_configurations])
    md_energies = np.array([md_output.spe[config - 1] for config in common_configurations])

    # print(qm_energies)
    # print(md_energies)
    # qm_energies = np.array(qm_energies)
    # md_energies = np.array(md_energies)

    # Create a DataFrame from the arrays
    df = pd.DataFrame({"QM Energies (kcal/mol)": qm_energies, "MD Energies (kcal/mol)": md_energies})

    # Write the DataFrame to an Excel and CSV files
    df.to_excel("energies.xlsx", index=True)
    df.to_csv("energies.csv", index=True)

    qm_energies = qm_energies - qm_energies.min()
    differences = qm_energies - md_energies

    average_difference = differences.mean()
    md_energies_shifted = md_energies + average_difference

    # Calculate the trendline (linear regression)
    coefficients = np.polyfit(md_energies_shifted, qm_energies, 1)
    polynomial = np.poly1d(coefficients)
    trendline = polynomial(md_energies_shifted)

    # Calculate R^2 value
    correlation_matrix = np.corrcoef(qm_energies, trendline)
    correlation_xy = correlation_matrix[0, 1]
    r_squared = correlation_xy**2

    # Create scatter plot
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=md_energies_shifted, y=qm_energies, mode="markers", name="Data"))
    fig.add_trace(go.Scatter(x=md_energies_shifted, y=trendline, mode="lines", name="Trendline"))
    fig.add_trace(
        go.Scatter(
            x=[min(md_energies_shifted), max(md_energies_shifted)],
            y=[min(qm_energies), max(qm_energies)],
            mode="lines",
            name="Linearity",
            line=dict(dash="dash"),
        )
    )

    # Adding the trendline equation and R^2 value to the plot
    trendline_eq_sign = "" if coefficients[1] < 0 else "+"
    fig.add_annotation(
        align="left",
        x=max(md_energies_shifted),
        y=min(qm_energies),
        text=f"Trendline: y(x)={coefficients[0]:.2f}x {trendline_eq_sign} {coefficients[1]:.2f} <br>R<sup>2</sup>={r_squared:.2f}",
        # text=f"Trendline: y={coefficients[0]:.2f}x+{coefficients[1]:.2f}\n$R^2$={r_squared:.2f}",
        showarrow=False,
        font=dict(size=18),
        xanchor="right",
    )

    # Additional figure settings
    fig.update_layout(
        # title="Scatter Plot with Trendline, Equation, and R²",
        xaxis_title="MD energies (kcal/mol)",
        yaxis_title="QM energies (kcal/mol)",
        font=dict(size=18),  # Set the font size here
    )

    # Save the figure
    file_path = "scatter_plot_with_trendline.png"  # Specify your desired file path and name here
    pio.write_image(fig, file_path)

    # print([(qm, md) for qm, md in zip(qm_energies.tolist(), md_energies_shifted.tolist())])
