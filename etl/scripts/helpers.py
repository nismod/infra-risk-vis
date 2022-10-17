"""
Standalone Snakemake helpers
"""
from typing import List

import pandas as pd


def hazard_layer_name_formatted(layer_name: str, hazard_layers: pd.DataFrame) -> str:
    """Generate the tileserver-formatted name string for the given layer"""
    try:
        print("getting hazard layer formatted name")
        name = (
            "{type}__rp_{rp}__rcp_{rcp}__epoch_{epoch}__conf_{confidence}.tif".format(
                type=hazard_layers[hazard_layers.key == layer_name].hazard.iloc[0],
                rp=hazard_layers[hazard_layers.key == layer_name].rp.iloc[0],
                rcp=hazard_layers[hazard_layers.key == layer_name].rcp.iloc[0],
                epoch=hazard_layers[hazard_layers.key == layer_name].epoch.iloc[0],
                confidence=hazard_layers[
                    hazard_layers.key == layer_name
                ].confidence.iloc[0],
            )
        )
        print("formatted name is", name)
        return name
    except IndexError as e:
        print(f"Could not find {layer_name} in hazard layers.")
        raise e
    finally:
        print("got hazard layer formatted name")


def generate_final_hazard_layer_names(hazard_layers: pd.DataFrame) -> List[str]:
    """Generate a list of all the expected final filenames"""
    names = []
    for layer_name in hazard_layers.key.to_list():
        names.append(hazard_layer_name_formatted(layer_name, hazard_layers))
    print("Generated output names:", names)
    return names
