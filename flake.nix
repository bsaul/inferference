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
    flake-utils.lib.eachDefaultSystem (system: 
    
    
    let

      package = "infeference";
      version = "1.0.3";
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

      # Run with
      # nix build .#cran
      packages.cran = pkgs.stdenv.mkDerivation {
          name = "cran";
          version = version;
          src = gitignoreSource ./.;
          buildInputs = [ 
              pkgs.R 
              pkgs.pandoc
              pkgs.rPackages.devtools
              pkgs.qpdf
              pkgs.texlive.combined.scheme-full
            ] ++ RDeps ;
          doCheck = true;
          buildPhase = ''
            ${pkgs.R}/bin/Rscript -e 'devtools::document()'
            ${pkgs.R}/bin/R CMD build .
          '';
          # NOTE: 
          # Not all checks will pass because R CMD check does stuff
          # This one gives warning:
          # * checking R files for syntax errors ... WARNING
          #  OS reports request to set locale to "en_US.UTF-8" cannot be honored
          # I don't know how to set the locale in a derivation 
          # (or if that's even possible)
          # Others are notes (e.g.):
          # * Found the following (possibly) invalid URLs:
          # R CMD check wants to access the interwebs but nix no like that.
          checkPhase = ''
            ${pkgs.R}/bin/R CMD check $(ls -t . | head -n1) --as-cran
          '';
          installPhase  = ''
            mkdir -p $out
            cp ${package}_${version}.tar.gz $out
            cp -r ${package}.Rcheck/ $out/logs
          '';
      };
      

      packages.default = pkgs.rPackages.buildRPackage {
        name = package;
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
