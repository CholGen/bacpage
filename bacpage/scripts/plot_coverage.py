import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Description"
    )

    # Initialize optional arguments
    parser.add_argument( "--input", type=str, help="Input depth file. " )
    parser.add_argument( "--output", type=str, help="output" )
    parser.add_argument( "--min-depth", type=int, help="Minimum required depth to call a position" )
    parser.add_argument( "--bin-size", type=int, help="Width of bins to " )

    args = parser.parse_args()

    depth = pd.read_csv( args.input, names=["chromosome", "position", "depth"], sep="\t" )
    depth["bin"] = (depth["position"] / args.bin_size).astype( int )
    depth = depth.groupby( ["chromosome", "bin"] )["depth"].agg( ["mean", "max", "min"] )
    depth = depth.reset_index()

    min_depth = args.min_depth

    chromosomes = depth["chromosome"].unique()

    fig, ax = plt.subplots( dpi=200, figsize=(10, 4 * len( chromosomes )), nrows=len( chromosomes ), sharey=True )
    ticks = mticker.FuncFormatter( lambda x, pos: '{0:,.0f}'.format( x * args.bin_size ) )

    leg = [
        Line2D( [0], [0], linestyle='-', marker=None, label=f"Average reads", markersize=0, linewidth=1 ),
        Line2D(
            [0], [0], linestyle='--', marker=None, color="black", label=f"Required depth ({min_depth})", markersize=0,
            linewidth=1
        )
    ]

    ax[0].legend( handles=leg, loc="upper left", handletextpad=1, frameon=False, fontsize=10, handlelength=1 )

    for i, name in enumerate( chromosomes ):
        data = depth.loc[depth["chromosome"] == name]
        ax[i].plot( "bin", "mean", data=data, zorder=10 )
        ax[i].fill_between( "bin", "min", "max", alpha=0.2, data=data, zorder=5 )
        ax[i].axhline( y=min_depth, color="black", linestyle="--", linewidth=1 )
        ax[i].xaxis.set_major_formatter( ticks )
        ax[i].set_xlabel( "Position (bp)" )
        ax[i].set_ylabel( "Reads" )
        ax[i].set_xlim( 0 )
        ax[i].set_ylim( 0 )
        ax[i].set_title( name )
        ax[i].grid( which="both", axis="both", linewidth=1, color="#F1F1F1", zorder=1 )

    plt.tight_layout()
    plt.savefig( args.output )
