{
  pkgs ? import (fetchTarball {
    url = "https://github.com/NixOS/nixpkgs/archive/refs/heads/release-23.05.tar.gz";
    sha256 = "sha256:0xhqjli4m9wkzv7xhs6fr1iajdjbv7xnj0bwvwldq9s6arlwkhj3";
  }) { },
  ...
}:
let
  python = pkgs.python39;
  src = pkgs.fetchurl {
    url = "https://files.pythonhosted.org/packages/63/ae/70fc292ccf056a9945e0f4591c29bbc5f018967e06eeb3f5a9449e0defff/catboost-1.2.5-cp39-cp39-macosx_11_0_universal2.whl";
    sha256 = "ed9521ce617f370dd7f8214f6b5ea2cdb2bfb47a2be0fb19a03b7c9d4399b4d1";
  };
in
python.pkgs.buildPythonPackage {
  pname = "catboost";
  version = "1.2.5";
  inherit src;
  format = "other";
  dontUnpack = true;
  installPhase = ''
    mkdir -p $out/lib/python3.9/site-packages
    ${pkgs.unzip}/bin/unzip -q $src -d $out/lib/python3.9/site-packages
  '';
  propagatedBuildInputs = with python.pkgs; [
    graphviz
    matplotlib
    numpy
    pandas
    plotly
    scipy
    six
  ];
}
