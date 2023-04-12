{
  inputs = {
    nixpkgs.url = "nixpkgs";
    flake-utils.url = "github:numtide/flake-utils";
  };
  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let pkgs = nixpkgs.legacyPackages.${system};
      in {
        packages = {
          codebook-gen = with pkgs;
            stdenv.mkDerivation rec {
              pname = "codebook-gen";
              version = "1.0";

              src = ./src-codebook;

              nativeBuildInputs = [ autoreconfHook gmp ];
            };
          sage = pkgs.sage.override {
            extraPythonPackages = ps: [
              (ps.buildPythonPackage rec {
                pname = "prtpy";
                version = "0.8.2";
                src = pkgs.fetchFromGitHub {
                  owner = "erelsgl";
                  repo = pname;
                  rev = "2a4e3c0077b0d71d30a04a1e273ed8986a28843e";
                  hash = "sha256-4Hxm0qejapEThYyl7PPxXi6MOQYACqhLkBRivQmOTt0=";
                };
                postPatch = ''
                  # Remove all version pinning (E.g., tornado==5.1.1 -> tornado)
                  sed -i 's/>=.*//g' requirements.txt
                '';
                doCheck = false;
                propagatedBuildInputs = with ps; [ numpy scipy mip ];
                meta.license = pkgs.lib.licenses.mit;
              })
            ];
            requireSageTests = false;
          };
        };

        apps = {
          sage = {
            type = "app";
            program = "${self.packages.${system}.sage}/bin/sage";
          };
        };
      });
}
