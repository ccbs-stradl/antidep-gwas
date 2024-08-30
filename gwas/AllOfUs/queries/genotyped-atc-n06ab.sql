SELECT
        d_exposure.person_id,
        d_exposure.drug_concept_id,
        d_standard_concept.concept_name as standard_concept_name,
        d_standard_concept.concept_code as standard_concept_code,
        d_standard_concept.vocabulary_id as standard_vocabulary,
        d_exposure.drug_exposure_start_datetime 
    FROM
        ( SELECT
            * 
        FROM
            `drug_exposure` d_exposure 
        WHERE
            (
                drug_concept_id IN (SELECT
                    DISTINCT ca.descendant_id 
                FROM
                    `cb_criteria_ancestor` ca 
                JOIN
                    (SELECT
                        DISTINCT c.concept_id       
                    FROM
                        `cb_criteria` c       
                    JOIN
                        (SELECT
                            CAST(cr.id as string) AS id             
                        FROM
                            `cb_criteria` cr             
                        WHERE
                            concept_id IN (21604709)             
                            AND full_text LIKE '%_rank1]%'       ) a 
                            ON (c.path LIKE CONCAT('%.', a.id, '.%') 
                            OR c.path LIKE CONCAT('%.', a.id) 
                            OR c.path LIKE CONCAT(a.id, '.%') 
                            OR c.path = a.id) 
                    WHERE
                        is_standard = 1 
                        AND is_selectable = 1) b 
                        ON (ca.ancestor_id = b.concept_id)))  
                    AND (d_exposure.PERSON_ID IN (SELECT
                        distinct person_id  
                FROM
                    `cb_search_person` cb_search_person  
                WHERE
                    cb_search_person.person_id IN (SELECT
                        person_id 
                    FROM
                        `cb_search_person` p 
                    WHERE
                        has_whole_genome_variant = 1 ) 
                    AND cb_search_person.person_id IN (SELECT
                        person_id 
                    FROM
                        `cb_search_person` p 
                    WHERE
                        has_ehr_data = 1 ) )
            )) d_exposure 
    LEFT JOIN
        `concept` d_standard_concept 
            ON d_exposure.drug_concept_id = d_standard_concept.concept_id;