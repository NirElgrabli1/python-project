"""
Population Analysis Script

This script analyzes population distribution and forecasts for Israel,
with a focus on Haredi population trends from 2009 to 2060.

Main functionalities:
- Load population data from Excel files
- Create visualizations of population distribution
- Generate forecast plots
- Analyze population trends

Required libraries: pandas, matplotlib, numpy, openpyxl
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import openpyxl
import warnings
import matplotlib

# Configure matplotlib parameters for consistent plotting
matplotlib.rcParams['font.family'] = 'Arial'
matplotlib.rcParams['axes.labelsize'] = 10
matplotlib.rcParams['axes.unicode_minus'] = False
matplotlib.rcParams['text.usetex'] = False
matplotlib.rcParams['axes.formatter.use_locale'] = True
plt.rcParams['figure.figsize'] = (15, 10)

# Suppress openpyxl warnings
warnings.filterwarnings('ignore', category=UserWarning, module='openpyxl')

# Define color palette for visualizations
COLORS = {
    'haredi': '#2E86C1',  # Blue color for Haredi population
    'non_haredi': '#F4D03F',  # Yellow color for non-Haredi population
    'arab': '#E74C3C'  # Red color for total population
}


def load_population_data():
    """
    Load population distribution and forecast data from Excel files.

    Attempts to read two Excel files:
    1. 2023 population distribution
    2. Population forecast from 2009 to 2060

    Returns:
    - DataFrame with 2023 population distribution
    - DataFrame with population forecast
    - Returns None if files cannot be loaded

    Expected Excel files:
    - '2023 פילוג אוכלוסיה.xlsx'
    - '2009 עד 2060 אוכלוסיה חרידית לעומת השאר.xlsx'
    """
    try:
        df_2023 = pd.read_excel('2023 פילוג אוכלוסיה.xlsx',
                                skiprows=2,
                                engine='openpyxl')
        df_2023 = df_2023[df_2023.iloc[:, 0].notna()]
        df_2023 = df_2023[~df_2023.iloc[:, 0].isin(['סך הכול', 'גברים', 'נשים'])]
        df_forecast = pd.read_excel('2009 עד 2060 אוכלוסיה חרידית לעומת השאר.xlsx',
                                    skiprows=2,
                                    engine='openpyxl')
        df_forecast = df_forecast[df_forecast.iloc[:, 1] == 'סך הכול']

        print("נתונים נטענו בהצלחה")
        return df_2023, df_forecast

    except FileNotFoundError:
        print("שגיאה: קובץ לא נמצא")
    except Exception as e:
        print(f"שגיאה בטעינת נתונים: {e}")

    return None, None


def create_age_pyramid(df):
    # Filter rows with age data
    mask = df.iloc[:, 0].str.contains(r'\d|\+', na=False)
    df_clean = df[mask].copy()
    df_clean = df_clean[~df_clean.iloc[:, 0].str.contains('סך הכול')]

    # Get unique age groups and sort them
    unique_ages = df_clean.iloc[:, 0].unique()
    age_order = sorted(unique_ages,
                       key=lambda x: int(x.split('-')[0].replace('+', '99')))  # Not reversed - youngest at bottom

    # Create aggregated data for each age group
    haredi_pct = []
    non_haredi_pct = []

    for age in age_order:
        age_data = df_clean[df_clean.iloc[:, 0] == age].iloc[0]
        haredi_pct.append(float(age_data.iloc[6]))
        non_haredi_pct.append(float(age_data.iloc[8]))

    y_pos = np.arange(len(age_order))

    # Create the plot
    fig, ax = plt.subplots(figsize=(15, 12))

    # Create the bars
    haredi_bars = ax.barh(y_pos, [-x for x in haredi_pct], color=COLORS['haredi'], label='םידרח', alpha=0.8)
    non_haredi_bars = ax.barh(y_pos, non_haredi_pct, color=COLORS['non_haredi'], label='םידרח אל', alpha=0.8)

    # Add percentage labels
    for i, v in enumerate(haredi_pct):
        if v >= 0.5:
            ax.text(-v - 0.3, i, f'{v:.1f}%', ha='right', va='center', fontsize=9)
    for i, v in enumerate(non_haredi_pct):
        if v >= 0.5:
            ax.text(v + 0.3, i, f'{v:.1f}%', ha='left', va='center', fontsize=9)

    # Customize the plot
    ax.set_yticks(y_pos)
    ax.set_yticklabels(age_order)
    ax.set_title('םידרח אל וא םידרח :םיאליג תוגלפתה, 2023', fontsize=14, pad=20)
    ax.set_xlabel('םיזוחא')

    # Add middle line
    ax.axvline(x=0, color='black', linewidth=0.5)

    # Add grid
    ax.grid(True, alpha=0.3, linestyle='--')

    # Add legend
    ax.legend(loc='upper right')

    # Set axis limits
    max_val = max(max(haredi_pct), max(non_haredi_pct))
    ax.set_xlim(-max_val - 5, max_val + 5)
    plt.tight_layout()
    plt.show()


def create_forecast_plot(df):
    # Create plot with white background
    fig, ax = plt.subplots(figsize=(15, 10), facecolor='white')

    # Prepare data
    years = df.iloc[:, 0].astype(int)
    total = df.iloc[:, 2].astype(float) / 1_000_000
    haredi = df.iloc[:, 6].astype(float) / 1_000_000
    non_haredi = df.iloc[:, 8].astype(float) / 1_000_000

    # Create the lines with markers
    ax.plot(years, haredi, '-o', color='#1f77b4', linewidth=2, label='םידרח', markersize=6)
    ax.plot(years, non_haredi, '-o', color='#ffcc00', linewidth=2, label='םידרח אל', markersize=6)
    ax.plot(years, total, '-o', color='#ff3333', linewidth=2, label='לכה ךס', markersize=6)

    # Add labels for last point only
    ax.text(years.iloc[-1], haredi.iloc[-1], f'{haredi.iloc[-1]:.1f}M', ha='left', va='bottom')
    ax.text(years.iloc[-1], non_haredi.iloc[-1], f'{non_haredi.iloc[-1]:.1f}M', ha='left', va='bottom')
    ax.text(years.iloc[-1], total.iloc[-1], f'{total.iloc[-1]:.1f}M', ha='left', va='bottom')

    # Set title and labels
    ax.set_title('0602-9002 תייסולכוא תיזחת', fontsize=14, pad=20)
    ax.set_xlabel('הנש', fontsize=12)
    ax.set_ylabel('(םינוילימב) הייסולכוא', fontsize=12)

    # Add grid
    ax.grid(True, alpha=0.3)

    # Set y axis range
    ax.set_ylim(0, max(total) * 1.1)

    # Create custom legend
    ax.legend(['םידרח םידוהי', 'םידרח אל םידוהי', 'לארשי תייסולכוא ללכ'],
              loc='upper left', fontsize=10)

    plt.tight_layout()
    plt.show()


def create_percentage_trend(df):
    """
    Create a bar plot showing the percentage trend of Haredi population.

    Visualization details:
    - Percentage of Haredi population over time
    - Gradient-colored bars representing percentage changes
    - Linear trend line to show overall population trend
    - Percentage labels on each bar

    Args:
    - df (pandas.DataFrame): Population forecast data
    """
    plt.figure(figsize=(15, 8), dpi=100)
    ax = plt.gca()

    # הכנת הנתונים
    years = df.iloc[:, 0].astype(int)
    total = df.iloc[:, 2].astype(float)
    haredi = df.iloc[:, 6].astype(float)
    percentages = (haredi / total) * 100

    # יצירת צבעים מדורגים לעמודות
    gradient_colors = plt.cm.Blues(np.linspace(0.5, 0.9, len(years)))

    # יצירת העמודות
    bars = ax.bar(years, percentages, color=gradient_colors,
                  width=2, alpha=0.7)

    # הוספת קו מגמה
    z = np.polyfit(range(len(years)), percentages, 1)
    p = np.poly1d(z)
    plt.plot(years, p(range(len(years))), '--',
             color='#E74C3C', linewidth=2,
             label='המגמ וק')

    # הוספת תוויות אחוזים מעל העמודות
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2., height + 0.3,
                f'{height:.1f}%',
                ha='center', va='bottom',
                fontsize=10,
                fontweight='bold')

    # עיצוב הגרף
    ax.grid(axis='y', linestyle='--', alpha=0.3)
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # כותרות וצירים
    plt.title('הייסולכואה ללכמ תידרחה הייסולכואה זוחא',
              pad=20, fontsize=14, fontweight='bold')
    plt.xlabel('הנש', fontsize=12)
    plt.ylabel('םיזוחא', fontsize=12)

    # הגדרת גבולות צירים
    ax.set_ylim(0, max(percentages) * 1.15)

    # הוספת מקרא
    plt.legend(loc='upper left')

    plt.tight_layout()
    plt.show()


def print_analysis(df):
    """
    Print a summary analysis of the population forecast.

    Calculates and displays key population statistics:
    - Final forecast year
    - Total population
    - Haredi population
    - Percentage of Haredi population

    Args:
    - df (pandas.DataFrame): Population forecast data
    """
    final_year = df.iloc[-1, 0]
    final_total = df.iloc[-1, 2] / 1_000_000
    final_haredi = df.iloc[-1, 6] / 1_000_000
    final_pct = (final_haredi / final_total) * 100

    print("\nניתוח האוכלוסייה החרדית בישראל")
    print("=" * 40)
    print(f"\nתחזית לשנת {final_year}")
    print(f"אוכלוסייה חרדית: {final_haredi:.2f} מיליון")
    print(f"אחוז מהאוכלוסייה: {final_pct:.1f}%")


def main():
    """
    Main function to orchestrate the population data analysis process.

    Workflow:
    1. Load population data from Excel files
    2. Create age pyramid visualization
    3. Generate population forecast plot
    4. Create percentage trend visualization
    5. Print population analysis summary
    """
    print("טוען קבצי נתונים...")
    pop_2023, forecast = load_population_data()

    if pop_2023 is not None and forecast is not None:
        print("\nיוצר תרשימים...")

        create_age_pyramid(pop_2023)
        create_forecast_plot(forecast)
        create_percentage_trend(forecast)

        print_analysis(forecast)

    else:
        print("\nלא ניתן להמשיך. אנא ודא שהקבצים הבאים קיימים באותו תיקייה:")
        print("1. '2023 פילוג אוכלוסיה.xlsx'")
        print("2. '2009 עד 2060 אוכלוסיה חרידית לעומת השאר.xlsx'")


if __name__ == "__main__":
    main()