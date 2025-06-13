{
  pkgs ? import (fetchTarball {
    url = "https://github.com/NixOS/nixpkgs/archive/refs/heads/release-23.05.tar.gz";
    sha256 = "sha256:0xhqjli4m9wkzv7xhs6fr1iajdjbv7xnj0bwvwldq9s6arlwkhj3";
  }) { },
  python ? pkgs.python39,
  ...
}:
python.pkgs.buildPythonPackage rec {
  pname = "logomaker";
  version = "0.8";
  src = pkgs.fetchPypi {
    inherit pname version;
    sha256 = "sha256-2MdQGn1teWHNaOWkTpOQAOvxsMQZegyRmDUeHWgdP20=";
  };
  pyproject = true;
  # dontUnpack = true;
  # installPhase = ''
  #   mkdir -p $out/lib/python3.9/site-packages
  #   ${pkgs.unzip}/bin/unzip -q $src -d $out/lib/python3.9/site-packages
  # '';
  propagatedBuildInputs = with python.pkgs; [
    pandas
    numpy
    matplotlib
  ];
}
