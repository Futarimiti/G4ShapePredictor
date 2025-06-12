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
        old-scikit-learn-overlay = self: super: {
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
          overlays = [ old-scikit-learn-overlay ];
        };
        python = pkgs.python39;
      in
      {
        # packages = rec {
        #   G4ShapePredictor = python.pkgs.buildPythonApplication {
        #     pname = "G4ShapePredictor";
        #     version = "1.0.0";
        #     pyproject = true;
        #     propagatedBuildInputs = with python.pkgs; [ setuptools ];
        #     src = ./.;
        #   };
        #   default = G4ShapePredictor;
        # };
        # apps = rec {
        #   {{cookiecutter.package_name}} = flake-utils.lib.mkApp {
        #     drv = self.packages.${system}.{{cookiecutter.package_name}};
        #   };
        #   default = {{cookiecutter.package_name}};
        # };
        devShells = {
          default = pkgs.mkShell {
            packages = [
              (python.withPackages (
                pypkgs: with pypkgs; [
                  pysimplegui
                  tkinter
                  biopython
                  pandas
                  numpy
                  scikit-learn
                  lightgbm
                  xgboost
                  # catboost # broken on darwin

                  ipython
                ]
              ))
            ];
          };
        };
      }
    );
}
