"""Download the public NHANES XPT components used in the analysis."""
from pathlib import Path
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
RAW.mkdir(parents=True, exist_ok=True)

FILES = {
    "2017": ["P_DEMO", "P_BMX", "P_BIOPRO", "P_CBC", "P_LUX", "P_BPXO",
             "P_GHB", "P_HDL", "P_TRIGLY", "P_DIQ", "P_BPQ", "P_ALQ", "P_HEQ", "P_HEPBD", "P_HEPC"],
    "2021": ["DEMO_L", "BMX_L", "BIOPRO_L", "CBC_L", "LUX_L", "BPXO_L",
             "GHB_L", "HDL_L", "TRIGLY_L", "DIQ_L", "BPQ_L", "ALQ_L", "HEQ_L", "HEPBD_L", "HEPC_L"],
}
BASE = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{cycle}/DataFiles/{name}.XPT"


def download(url: str, destination: Path) -> None:
    request = Request(url, headers={"User-Agent": "NHANES-reproducibility-script/1.0"})
    with urlopen(request, timeout=120) as response, destination.open("wb") as output:
        output.write(response.read())
    if destination.stat().st_size < 1000:
        destination.unlink(missing_ok=True)
        raise RuntimeError(f"Downloaded file is unexpectedly small: {url}")


def main() -> None:
    for cycle, names in FILES.items():
        for name in names:
            destination = RAW / f"{name}.xpt"
            if destination.exists() and destination.stat().st_size >= 1000:
                print(f"exists  {destination.name}")
                continue
            url = BASE.format(cycle=cycle, name=name)
            print(f"download {url}")
            download(url, destination)


if __name__ == "__main__":
    main()
