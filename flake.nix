{
  description = "inferference";
  nixConfig = {
    bash-prompt = "> ";
  };
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-24.11";
    flake-utils.url = "github:numtide/flake-utils";
    gitignore = {
      url = "github:hercules-ci/gitignore.nix";
      # Use the same nixpkgs
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  outputs = { self, nixpkgs, flake-utils, gitignore }:
    flake-utils.lib.eachDefaultSystem (system: let
      pkgs = nixpkgs.legacyPackages.${system};

      RDeps = with pkgs.rPackages; 
        [ numDeriv 
          lme4
          Formula 
          spaMM
          testthat
          knitr
          markdown
        ];

      inherit (gitignore.lib) gitignoreSource;

    in {

      packages.default = pkgs.rPackages.buildRPackage {
        name = "inferference";
        src  = gitignoreSource ./.;
        propagatedBuildInputs = RDeps;
      };

      devShells.default =  pkgs.mkShell {
        buildInputs = [

          # R and packages
          pkgs.R
          pkgs.rPackages.devtools
          pkgs.rPackages.languageserver

          # Documentation/writing tools
          pkgs.pandoc
        ] ++ RDeps;
      }; 
    });
}
