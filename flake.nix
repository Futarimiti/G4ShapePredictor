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
                  inherit (old) pname;
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
            # catboost # FIXME no support on darwin
          ];
      in
      {
        packages = rec {
          g4sp =
            let
              pythonEnv = python.withPackages (pypkgs: dependenciesFrom pypkgs);
            in
            pkgs.stdenv.mkDerivation {
              pname = "g4sp";
              version = "1.0.0";
              src = ./.;
              buildInputs = [ pythonEnv ];
              installPhase = ''
                mkdir -p $out/bin $out/share/g4sp
                cp -r "g4sp application code" $out/share/g4sp
                cat > $out/bin/g4sp << EOF
                #!${pkgs.bash}/bin/bash
                cd $out/share/g4sp
                exec ${pythonEnv}/bin/python "g4sp application code/G4ShapePredictor.py"
                EOF
                chmod +x $out/bin/g4sp
              '';
            };
          default = g4sp;
        };
        apps = rec {
          g4sp = flake-utils.lib.mkApp { drv = self.packages.${system}.g4sp; };
          default = g4sp;
        };
        devShells.default = pkgs.mkShell {
          packages = [ (python.withPackages (pypkgs: dependenciesFrom pypkgs ++ [ pypkgs.ipython ])) ];
        };
      }
    );
}
