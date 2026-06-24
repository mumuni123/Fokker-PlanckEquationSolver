from plot_beam_density_profile import main as beam_main
from plot_background_density_evolution import main as main_background

def main() -> None:
    main_background()
    beam_main()

if __name__ == "__main__":
    main()