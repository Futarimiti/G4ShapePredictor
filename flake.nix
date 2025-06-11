{
  description = "Accurate folding topologies predictor for DNA G4 in potassium (K+) buffer, based on putative quadruplex sequences";
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs";
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
        pkgs = nixpkgs.legacyPackages.${system};
        pythonVersion = "311";
        python = pkgs."python${pythonVersion}";
      in
      {
        # packages = rec {
        #   {{cookiecutter.package_name}} = python.pkgs.buildPython{{cookiecutter.app_or_pkg.capitalize()}} {
        #     pname = "{{cookiecutter.package_name}}";
        #     version = "{{cookiecutter.package_version}}";
        #     pyproject = true;
        #     propagatedBuildInputs = with python.pkgs; [ setuptools ];
        #     src = ./.;
        #   };
        #   default = {{cookiecutter.package_name}};
        # };
        # apps = rec {
        #   {{cookiecutter.package_name}} = flake-utils.lib.mkApp {
        #     drv = self.packages.${system}.{{cookiecutter.package_name}};
        #   };
        #   default = {{cookiecutter.package_name}};
        # };
        devShells = {
          default = pkgs.mkShell {
            packages =
              let
                pysimplegui =
                  with python.pkgs;
                  buildPythonPackage {
                    pname = "PySimpleGUI";
                    version = "4.60.5";
                    src = fetchTarball {
                      url = "https://deleted-packages.pypimirror.stablebuild.com/pysimplegui/PySimpleGUI-4.60.5.tar.gz?expires=1749609803478&signature=1d086905c6dabbd415dec36054f1c41d23f4f56b076c8d797a73fbaa643ce193";
                      sha256 = "1jpgighmkc3pw5n0g5i6cdjyczj5yp6h0yrmiwjbm3myhsf3rlrf";
                    };
                    doCheck = false;
                  };
              in
              [
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
                  ]
                ))
              ];
          };
        };
      }
    );
}
