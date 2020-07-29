
cd /content/drive/My Drive

pip install xlsxwriter

import xml.etree.ElementTree as ET
from itertools import combinations
import xlsxwriter

'''
XML Schema is like - 
Root : root.tag = {http://www.drugbank.ca}drugbank ; root.attrib = {exported on etc etc}
  Child : Drug
    ChildOfDrug: Id
    ChildOfDrug: Name
    .
    .
    .
    .
    .
    ChildOfDrug: 
  Child : Drug
'''

filename = "full_database.xml"
tree = ET.parse(filename)
root = tree.getroot()
root_str = root.tag.split('}')[0]+'}'

def get_numOfDrugs(element=root):
  return len(element)

def get_numOfChildren(element):
  return len(element)

def get_childByDrugName(name):
  found = False
  for child in root:
    if(name==child.findtext(root_str+'name')):
      found = True
      break
  return [found,child]

def get_listOfDrugNames(element=root): 
  return [child.findtext(root_str+'name') for child in element]

def get_listOfDrugbankIds(element=root):
  return [child[0].text for child in element]

def get_drugId(drugname):
  lst = get_childByDrugName(drugname)
  if (lst[0]):
    child = lst[1]
    return child[0].text
  else:
    lst = get_childByDrugName(drugname.split()[0])
    if (lst[0]):
      return lst[1].text
    else:
      print ("Drug Not Found")
      return "" 

def get_drugName(id):
  found = False
  for child in root:
    if(child[0].text == id):
      found = True
      break
  name = child.findtext(root_str+'name') if (found) else  ""
  return name

def get_attribOfDrug(name,attrib):
  #Valid Attribs = {'name','description','cas-number','unii','state' etc}
  #There is no check for attrib - enter valid attribs
  #Incase an empty string is returned means the attrib is incorrect
  lst = get_childByDrugName(name)
  if (lst[0]):
    return lst[1].findtext(root_str+attrib)
  else:
    print ("Drug Not Found")
    return ""

def get_drugGroup(name):
  lst = get_childByDrugName(name) 
  if (lst[0]):
    return [c.text for c in lst[1].find(root_str+'groups')]
  else:
    print ("Drug Not Found")
    return ""

def get_listOfApprovedDrugs():
  drug_list = []
  for child in root:
    group = get_drugGroup(child.findtext(root_str+'name'))
    if ('approved' in group):
      drug_list.append(child.findtext(root_str+'name'))
  return drug_list

def get_drugInteractions(drugname):
  interaction_list = []
  lst = get_childByDrugName(drugname)
  if(lst[0]):
    element = lst[1].find(root_str+'drug-interactions')
    return [intr[1].text for intr in element] #index 0 - gives the ids of the interacting drugs
  else:
    print ("DrugName is Incorrect")
    return ""

def check_drugInteraction(drug1,drug2):
  c1 = get_drugId(drug1)
  c2 = get_drugId(drug2)
  if(c1 == "" or c2 == ""):
    print ("Check the Drug Names again")
    return ""
  else:
    flag = True if (drug2 in get_drugInteractions(drug1)) else False
    return flag

def get_SMPDBProteins(drugname): 
  proteins = []
  lst = get_childByDrugName(drugname)
  if (lst[0]):
    elem = lst[1].find(root_str+'pathways')
    if (len(elem)!=0):
      elem = elem[0].find(root_str+'enzymes')
      proteins = [elem[i].text for i in range(len(elem))]
    else:
      print ("No SMPDB Target Interactions listed")
  else:
    print ("Check the Drug Name again")
  return proteins

def get_Proteins(drugname,typeOfProtein): #Type = 'transporters' Or 'enzymes' Or 'targets'
  proteins=[]
  c = get_childByDrugName(drugname)
  if (c[0]):
    c = c[1].find(root_str+typeOfProtein)
    if (len(c)>0):
      for child in c:
        p = child.find(root_str+'polypeptide')
        if (not (p is None)):
          proteins.append(p.attrib['id'])
    else:
      print ("%s : %s :No Listed Porteins"%(drugname,typeOfProtein))
  else:
    print ("Recheck the DrugName")
  return proteins

def check_typeOfInteraction(drug1,drug2,typeOfInteraction): #Type = 'transporters' Or 'enzymes' Or 'targets'
  c1 = get_childByDrugName(drug1)
  if (c1[0]):
    c2 = get_childByDrugName(drug2)
    if (c2[0]):
      if (check_drugInteraction(drug1,drug2)):
        p1 = get_Proteins(drug1,typeOfInteraction)
        p2 = get_Proteins(drug2,typeOfInteraction)
        intr = set(p1).intersection(set(p2))
        flag = True if (len(intr)>0) else False
        return flag
      else :
        print ("Drugs do not interact")
        return ""
    else:
      print ("Drug2 not found")
      return ""
  else :
    print ("Drug1 not Found")
    return ""

def get_foodInteractions(drugname):
  food_interactions = []
  c = get_childByDrugName(drugname)
  if (c[0]):
    c = c[1].find(root_str+'food-interactions')
    if(len(c)>0):
      food_interactions = [child.text for child in c]
  else:
    print("Drug not found")
  return food_interactions

def get_drugLinks(druglist):
  links = []
  drugs_done=[]
  for elem in druglist:
    intr = get_drugInteractions(elem)
    approved_intr = list(set(intr).intersection(set(druglist)))
    for d in drugs_done:
      if (d in approved_intr):
        approved_intr.remove(d)
    drugs_done.append(elem)
    for d in approved_intr:
      links.append([elem,d])
  return links

def get_ChEMBL_id(drugname):
  id = ""
  c = get_childByDrugName(drugname)
  if (c[0]):
    c = c[1].find(root_str+'external-identifiers')
    for child in c:
      if (child[0].text == "ChEMBL"):
        break
    id = child[1].text
  else :
    print ("Recheck the drugname")
  return id

def generate_drugbank_ChEMBL_id_mapping_file():
  workbook = xlsxwriter.Workbook('Drugbank-ChEMBL.xlsx')
  worksheet = workbook.add_worksheet()
  row,col=0,0
  for child in root:
    name = child.find(root_str+'name').text
    d_id = child.find(root_str+'drugbank-id').text
    c = child.find(root_str+'external-identifiers')
    for elem in c:
      if (elem[0].text == "ChEMBL"):
        c_id = elem[1].text
        break
    else:
      c_id = "--"

    worksheet.write(row, col, name)
    worksheet.write(row, col+1, d_id)
    worksheet.write(row, col+2, c_id)
    row += 1
  workbook.close()

def generate_drugbank_CHEBI_id_mapping_file():
  workbook = xlsxwriter.Workbook('Drugbank-CHEBI.xlsx')
  worksheet = workbook.add_worksheet()
  row,col=0,0
  for child in root:
    name = child.find(root_str+'name').text
    d_id = child.find(root_str+'drugbank-id').text
    c = child.find(root_str+'external-identifiers')
    for elem in c:
      if (elem[0].text == "ChEBI"):
        c_id = elem[1].text
        break
    else:
      c_id = "--"

    worksheet.write(row, col, name)
    worksheet.write(row, col+1, d_id)
    worksheet.write(row, col+2, c_id)
    row += 1
  workbook.close()

#We can use the APIs this way -- 
generate_drugbank_CHEBI_id_mapping_file()
c = get_childByDrugName('Stanolone')
id = get_ChEMBL_id('Dexamethasone')

#Gives the interactions between the drugs mentioned in a list
get_drugLinks(['Lepirudin','Acetaminophen','Ephedra sinica root','Dapagliflozin','Ibuprofen','Denileukin diftitox'])

get_drugGroup('Denileukin diftitox')

get_drugInteractions('Ibuprofen')

check_typeOfInteraction('Phenelzine','Acetaminophen','enzymes')

approved_drugs = get_listOfApprovedDrugs()

get_numOfDrugs()

get_drugId('Hydrochlorothiazide')