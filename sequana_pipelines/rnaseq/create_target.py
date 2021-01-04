import glob
import click

version ="1.0"
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=version)
@click.option('--feature-counts-directory', default="feature_counts")
def main(feature_counts_directory):
    """Create target file for the RNADiff analysis"""
    filenames = glob.glob("{}/*_feature.out".format(
            feature_counts_directory))
    print("label\tfiles\tcondition\treplicat")
    for filename in filenames:
        label = filename.split("/")[-1].replace("_feature.out", "")
        filename = filename.split("/")[-1]
        print("{}\t{}\t{}\t{}".format(label, filename, label, "X"))


if __name__ == "__main__": #pragma: no cover
    main()

