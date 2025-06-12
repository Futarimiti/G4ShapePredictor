{
  description = "Accurate folding topologies predictor for DNA G4 in potassium (K+) buffer, based on putative quadruplex sequences";
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/release-23.05"; # python39 support removed onwards
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
      ...
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        python-packages-overlay = self: super: {
          pythonPackagesExtensions = super.pythonPackagesExtensions ++ [
            (pyself: pysuper: {
              scikit-learn = pysuper.scikit-learn.overridePythonAttrs (old: rec {
                version = "1.0.2";
                src = pysuper.fetchPypi {
                  pname = old.pname;
                  inherit version;
                  sha256 = "sha256-tYcJWaVIS2FPJtMcpMF1JLGwMXUiGZ3JhcO0JW4DB2c=";
                };
              });
            })
          ];
        };
        pkgs = import nixpkgs {
          inherit system;
          overlays = [ python-packages-overlay ];
        };
        python = pkgs.python39;
        dependenciesFrom =
          pypkgs: with pypkgs; [
            pysimplegui
            tkinter
            biopython
            pandas
            numpy
            scikit-learn
            lightgbm
            xgboost
            # catboost # no support on darwin
          ];
      in
      {
        packages = rec {
          G4ShapePredictor = throw "TODO: packages.G4ShapePredictor";
          default = G4ShapePredictor;
        };
        apps = rec {
          G4ShapePredictor = flake-utils.lib.mkApp { drv = self.packages.${system}.G4ShapePredictor; };
          default = G4ShapePredictor;
        };
        devShells.default = pkgs.mkShell {
          packages = [ (python.withPackages (pypkgs: dependenciesFrom pypkgs ++ [ pypkgs.ipython ])) ];
        };
      }
    );
}
