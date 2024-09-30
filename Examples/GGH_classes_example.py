from GGH_crypto import GGHCryptosystem, GGHHNFCryptosystem

dimension = 100
GGH = GGHCryptosystem(dimension = dimension, debug=True)
GGH.encrypt() #You first need to encrypt

message = GGH.message

decrypted_message = GGH.decrypt() #and then decrypt

if decrypted_message == message:
    print("GGH decryption succeeded")

GGH_HNF = GGHHNFCryptosystem(dimension = dimension, debug=True)
GGH_HNF.encrypt() #You first need to encrypt

message = GGH_HNF.error #The message is stored in the error vector in GGH-HNF

decrypted_message = GGH_HNF.decrypt() #and then decrypt

if decrypted_message == message:
    print("GGH-HNF decryption succeeded")
else:
    print("Decryption failed, try with GGH_private = True or lower the alpha value")