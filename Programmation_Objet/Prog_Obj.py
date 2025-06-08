class Personne:
    compteur = 0

    def __init__(self, nom, age):
        self.nom = nom 
        self.age = age
        Personne.compteur += 1

    def se_presenter(self):
        print(f"Je m'appelle {self.nom} et j'ai {self.age} ans")

    @classmethod
    def afficher_compteur(cls):
        print(f"Nombre de personnes créées : {Personne.compteur}")


p = Personne("Alice", 30)
p.se_presenter()
p.afficher_compteur()

class Employe(Personne):
    def __init__(self, nom, age, salaire):
        super().__init__(nom, age)
        self.salaire = salaire

    def se_presenter(self):
        super().se_presenter()
        print(f"Je gagne {self.salaire} euros")

e = Employe("Charlie", 40, 3500)
e.se_presenter()

        

    