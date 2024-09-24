import pandas as pd
import requests
import os

def url_from_key(key):
    """
    Lookup an Aqueduct TIFF URL from our layers file by KEY value.
    """
    # Read the CSV file containing the layer information
    df = pd.read_csv("pipelines/aqueduct/layers.csv")

    # Find the row where the filename matches the given key
    layer = df[df['filename'] == f"{key}.tif"].squeeze()

    # Extract the URL from the row
    url = layer['url']
    return url

def download_file(url, output_path):
    """
    Download a file from the given URL and save it to the specified output path.
    """
    # Ensure the output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Make the HTTP request to download the file
    response = requests.get(url, stream=True)

    # Write the content to the output path
    with open(output_path, 'wb') as f:
        f.write(response.content)

    print(f"File saved to {output_path}")

def main():
    # Specify the key for which we want to download the file
    key = 'inuncoast_historical_wtsub_hist_rp0001_5'  # Replace this with the actual key you are testing

    # Get the URL from the key
    url = url_from_key(key)
    print(f"URL for {key}: {url}")

    # Define the output path where the file should be saved
    output_path = f"raster/raw/aqueduct/{key}.tif"

    # Download the file from the URL and save it to the output path
    download_file(url, output_path)

if __name__ == "__main__":
    main()
